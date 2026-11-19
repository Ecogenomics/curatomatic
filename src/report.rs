//! The review table for curators: one row per taxon that curatomatic changed or flagged.
//!
//! The report is computed from data already in memory and never changes the tree.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::fs;
use std::path::Path;

use strum::{Display, IntoEnumIterator};

use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::metadata::Metadata;
use crate::newick::{NewickTree, NodeId};
use crate::rank::Rank;
use crate::red::{Red, RedValue};
use crate::taxon::Taxon;
use crate::tree::{Move, Tree};

/// A move that multiplies the number of genomes by this much or more is a `large_gain`.
const LARGE_GAIN_FACTOR: usize = 5;

/// Nodes below this `q` are low for their rank (see [`Thresholds::q`]).
const LOW_Q: f64 = 0.5;

/// Separates multiple values in one column.
const LIST_SEPARATOR: &str = ",";

const HEADER: [&str; 13] = [
    "taxon",
    "rank",
    "action",
    "flag",
    "n_genomes_before",
    "n_genomes_after",
    "bootstrap_before",
    "bootstrap_after",
    "red_before",
    "red_after",
    "prior_taxon",
    "nested_in",
    "sister_taxa",
];

/// What curatomatic did to a taxon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Display)]
#[strum(serialize_all = "snake_case")]
pub(crate) enum Action {
    MovedUp,
    MovedDown,
    Novel,
}

/// Something about a taxon that a curator should review.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Display)]
#[strum(serialize_all = "snake_case")]
pub(crate) enum Flag {
    /// The taxon sits inside another taxon of the same rank.
    NestedInSameRank,
    /// A novel taxon with two or more genomes on a node with `q < 0.5`.
    NovelLowRed,
    /// The taxon joins its closest sister taxon of the same rank above `child_thresh`.
    SisterJoinPastChildThresh,
    /// The taxon is on a node with a bootstrap below `--min-bs`.
    UnsupportedLabel,
    /// A move that multiplied the number of genomes by 5 or more.
    LargeGain,
    /// An unmoved taxon on a supported node with `q < 0.5`, whose genomes from the previous
    /// release are all under one child, while the other children hold only new genomes.
    SplitCandidate,
    /// A novel taxon whose genomes had a taxon of this rank in the previous release.
    DisplacedFromPrior,
    /// A taxon on the input tree that is in no `gtdb_taxonomy` of the metadata.
    NewInDecoration,
}

/// The node that holds a taxon, as reported.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NodeStats {
    pub(crate) n_genomes: usize,
    pub(crate) bootstrap: Option<BootstrapValue>,
    pub(crate) red: RedValue,
}

/// One row of the report.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ReportRow {
    pub(crate) taxon: Taxon,
    pub(crate) rank: Rank,
    pub(crate) action: Option<Action>,
    pub(crate) flags: BTreeSet<Flag>,
    /// The node of the taxon on the input tree, a novel taxon has none.
    pub(crate) before: Option<NodeStats>,
    pub(crate) after: NodeStats,
    /// The previous-release taxa of this rank held by the genomes of a novel taxon.
    pub(crate) prior_taxa: BTreeSet<Taxon>,
    /// The closest taxon of the same rank that this taxon is nested inside.
    pub(crate) nested_in: Option<Taxon>,
    /// The closest taxa of the same rank that join this one past `child_thresh`.
    pub(crate) sister_taxa: BTreeSet<Taxon>,
}

/// Everything the report is computed from.
pub(crate) struct ReportInput<'a> {
    /// The curated tree, reloaded and still holding its RED values.
    pub(crate) tree: &'a Tree,
    /// The node of every taxon on the input tree.
    pub(crate) taxa_before: &'a BTreeMap<Taxon, NodeId>,
    /// The moves made by refinement.
    pub(crate) moves: &'a [Move],
    /// The novel taxa that were created.
    pub(crate) novel: &'a BTreeSet<Taxon>,
    /// The metadata, subset to the genomes on the tree.
    pub(crate) metadata: &'a Metadata,
    /// Every taxon in a `gtdb_taxonomy` of the full metadata file.
    pub(crate) prior_taxa: &'a BTreeSet<Taxon>,
    pub(crate) min_bs: BootstrapValue,
    pub(crate) red_dict: &'a Red,
}

/// The RED thresholds of a rank `r`.
#[derive(Debug, Clone, Copy)]
struct Thresholds {
    median: f64,
    /// The midpoint of the median of the parent rank and `median`.
    parent: f64,
    /// The midpoint of `median` and the median of the child rank.
    child: f64,
}

impl Thresholds {
    fn new(rank: Rank, red_dict: &Red) -> CuratomaticResult<Self> {
        let median = red_dict.get(rank)?;
        let parent = rank
            .parent()
            .ok_or(CuratomaticError::NoParentRank(rank.prefix()))?;
        let child = rank.child().ok_or_else(|| {
            CuratomaticError::Invariant(format!("the {}__ rank has no child rank", rank.prefix()))
        })?;
        Ok(Thresholds {
            median: median.value(),
            parent: red_dict.get(parent)?.midpoint(median).value(),
            child: red_dict.get(child)?.midpoint(median).value(),
        })
    }

    /// Where a RED sits within the rank: 0 at the parent threshold, 1 at the median.
    fn q(self, red: RedValue) -> f64 {
        (red.value() - self.parent) / (self.median - self.parent)
    }
}

/// The ranks that are refined, and flagged.
fn is_flagged_rank(rank: Rank) -> bool {
    matches!(
        rank,
        Rank::Phylum | Rank::Class | Rank::Order | Rank::Family | Rank::Genus
    )
}

/// Data derived from the structure and labels of the curated tree.
struct TreeIndex<'a> {
    tree: &'a Tree,
    /// The number of genomes (leaves that are not zombies) below each node.
    n_genomes: Vec<usize>,
}

impl<'a> TreeIndex<'a> {
    fn new(tree: &'a Tree) -> CuratomaticResult<Self> {
        let newick = tree.newick();
        let mut n_genomes = Vec::new();
        for node in newick.ids() {
            let n = if newick.is_leaf(node) {
                usize::from(
                    tree.label(node)?
                        .gid
                        .as_ref()
                        .is_some_and(|gid| !gid.is_zombie()),
                )
            } else {
                newick
                    .children(node)
                    .iter()
                    .map(|child| n_genomes.get(child.index()).copied().unwrap_or_default())
                    .sum()
            };
            n_genomes.push(n);
        }
        Ok(TreeIndex { tree, n_genomes })
    }

    fn newick(&self) -> &'a NewickTree {
        self.tree.newick()
    }

    fn n_genomes(&self, node: NodeId) -> usize {
        self.n_genomes
            .get(node.index())
            .copied()
            .unwrap_or_default()
    }

    fn stats(&self, node: NodeId) -> CuratomaticResult<NodeStats> {
        let label = self.tree.label(node)?;
        Ok(NodeStats {
            n_genomes: self.n_genomes(node),
            bootstrap: label.bootstrap,
            red: label.red,
        })
    }

    /// The genomes (not zombies) at or below `node`, which may be a leaf itself.
    fn genomes_under(&self, node: NodeId) -> CuratomaticResult<Vec<&'a GenomeId>> {
        let mut out = Vec::new();
        for desc in self.newick().descendants(node).chain([node]) {
            if let Some(gid) = &self.tree.label(desc)?.gid
                && !gid.is_zombie()
            {
                out.push(gid);
            }
        }
        Ok(out)
    }

    /// The named taxa of `rank` on the label of `node`.
    fn named(&self, node: NodeId, rank: Rank) -> CuratomaticResult<Vec<&'a Taxon>> {
        Ok(self
            .tree
            .label(node)?
            .taxonomy
            .as_ref()
            .map_or_default(|taxonomy| {
                taxonomy
                    .get_rank(rank)
                    .iter()
                    .filter(|t| t.is_named())
                    .collect()
            }))
    }

    fn is_ancestor(&self, ancestor: NodeId, node: NodeId) -> bool {
        let mut current = self.newick().parent(node);
        while let Some(cur) = current {
            if cur == ancestor {
                return true;
            }
            current = self.newick().parent(cur);
        }
        false
    }
}

/// The number of named taxa of a rank on each node, and in the subtree of each node.
struct RankCounts {
    own: Vec<usize>,
    subtree: Vec<usize>,
}

impl RankCounts {
    fn new(index: &TreeIndex<'_>, rank: Rank) -> CuratomaticResult<Self> {
        let newick = index.newick();
        let (mut own, mut subtree) = (Vec::new(), Vec::new());
        for node in newick.ids() {
            let n_own = index.named(node, rank)?.len();
            let n_below: usize = newick
                .children(node)
                .iter()
                .map(|child| subtree.get(child.index()).copied().unwrap_or_default())
                .sum();
            own.push(n_own);
            subtree.push(n_own.saturating_add(n_below));
        }
        Ok(RankCounts { own, subtree })
    }

    fn own(&self, node: NodeId) -> usize {
        self.own.get(node.index()).copied().unwrap_or_default()
    }

    fn subtree(&self, node: NodeId) -> usize {
        self.subtree.get(node.index()).copied().unwrap_or_default()
    }
}

/// The closest taxa of `rank` that join the taxon on `node`, if they join above `child_thresh`.
///
/// Walks up from `node` to the first ancestor with a taxon of the rank under another child, and
/// returns the highest taxa of the rank under those children. A nested taxon (one that reaches an
/// ancestor holding a taxon of the rank first) has no sister.
fn sister_join(
    index: &TreeIndex<'_>,
    counts: &RankCounts,
    rank: Rank,
    node: NodeId,
    thresholds: Thresholds,
) -> CuratomaticResult<BTreeSet<Taxon>> {
    let newick = index.newick();
    let mut out = BTreeSet::new();
    let mut child = node;
    let mut current = newick.parent(node);
    while let Some(ancestor) = current {
        if counts.own(ancestor) > 0 {
            break;
        }
        let other = counts
            .subtree(ancestor)
            .saturating_sub(counts.subtree(child));
        if other > 0 {
            if index.tree.label(ancestor)?.red.value() > thresholds.child {
                let mut stack: Vec<NodeId> = newick
                    .children(ancestor)
                    .iter()
                    .copied()
                    .filter(|&c| c != child)
                    .collect();
                while let Some(cur) = stack.pop() {
                    if counts.subtree(cur) == 0 {
                        continue;
                    }
                    let named = index.named(cur, rank)?;
                    if named.is_empty() {
                        stack.extend(newick.children(cur));
                    } else {
                        out.extend(named.into_iter().cloned());
                    }
                }
            }
            break;
        }
        child = ancestor;
        current = newick.parent(ancestor);
    }
    Ok(out)
}

/// True if the genomes below `node` that the metadata puts in `taxon` are all under one child,
/// and every other child holds only genomes that are new (have no prior taxonomy).
fn is_split_candidate(
    index: &TreeIndex<'_>,
    metadata: &Metadata,
    rank: Rank,
    taxon: &Taxon,
    node: NodeId,
) -> CuratomaticResult<bool> {
    let children = index.newick().children(node);
    if children.len() < 2 {
        return Ok(false);
    }
    let mut n_with_members = 0_usize;
    let mut others_new = true;
    for &child in children {
        let genomes = index.genomes_under(child)?;
        let has_members = genomes
            .iter()
            .any(|gid| metadata.taxa_of_rank(gid, rank).any(|t| t == taxon));
        if has_members {
            n_with_members = n_with_members.saturating_add(1);
        } else if genomes.is_empty() || genomes.iter().any(|gid| metadata.has_taxonomy(gid)) {
            others_new = false;
        }
    }
    Ok(n_with_members == 1 && others_new)
}

/// Build the rows of the report, sorted by rank and then taxon.
pub(crate) fn build_report(input: &ReportInput<'_>) -> CuratomaticResult<Vec<ReportRow>> {
    let index = TreeIndex::new(input.tree)?;
    let taxa_after = input.tree.taxon_nodes();
    let has_prior = !input.prior_taxa.is_empty();

    let moves: BTreeMap<&Taxon, &Move> = input.moves.iter().map(|m| (&m.taxon, m)).collect();

    let mut thresholds = BTreeMap::new();
    let mut counts = BTreeMap::new();
    for rank in Rank::iter().filter(|&rank| is_flagged_rank(rank)) {
        thresholds.insert(rank, Thresholds::new(rank, input.red_dict)?);
        counts.insert(rank, RankCounts::new(&index, rank)?);
    }

    let mut rows = Vec::new();
    for (taxon, &node) in &taxa_after {
        let rank = taxon.rank()?;
        let is_novel = input.novel.contains(taxon);
        let after = index.stats(node)?;

        let (action, before) = if is_novel {
            (Some(Action::Novel), None)
        } else if let Some(m) = moves.get(taxon) {
            let action = if index.is_ancestor(m.to, m.from) {
                Action::MovedUp
            } else {
                Action::MovedDown
            };
            (Some(action), Some(index.stats(m.from)?))
        } else {
            let from = input.taxa_before.get(taxon).copied().unwrap_or(node);
            (None, Some(index.stats(from)?))
        };

        let mut row = ReportRow {
            taxon: taxon.clone(),
            rank,
            action,
            flags: BTreeSet::new(),
            before,
            after,
            prior_taxa: BTreeSet::new(),
            nested_in: None,
            sister_taxa: BTreeSet::new(),
        };

        if taxon.is_named()
            && let (Some(&thresh), Some(rank_counts)) = (thresholds.get(&rank), counts.get(&rank))
        {
            let q = thresh.q(after.red);

            if let Some(outer) = input.tree.nested_inside(node, rank, taxon) {
                row.flags.insert(Flag::NestedInSameRank);
                row.nested_in = Some(outer.clone());
            }

            if is_novel && after.n_genomes >= 2 && q < LOW_Q {
                row.flags.insert(Flag::NovelLowRed);
            }

            row.sister_taxa = sister_join(&index, rank_counts, rank, node, thresh)?;
            if !row.sister_taxa.is_empty() {
                row.flags.insert(Flag::SisterJoinPastChildThresh);
            }

            if after.bootstrap.is_some_and(|bs| bs < input.min_bs) {
                row.flags.insert(Flag::UnsupportedLabel);
            }

            if let (Some(Action::MovedUp | Action::MovedDown), Some(before)) = (action, before)
                && before.n_genomes > 0
                && after.n_genomes >= before.n_genomes.saturating_mul(LARGE_GAIN_FACTOR)
            {
                row.flags.insert(Flag::LargeGain);
            }

            if action.is_none()
                && after.bootstrap.is_some_and(|bs| bs >= input.min_bs)
                && q < LOW_Q
                && is_split_candidate(&index, input.metadata, rank, taxon, node)?
            {
                row.flags.insert(Flag::SplitCandidate);
            }

            if is_novel {
                for gid in index.genomes_under(node)? {
                    row.prior_taxa
                        .extend(input.metadata.taxa_of_rank(gid, rank).cloned());
                }
                if !row.prior_taxa.is_empty() {
                    row.flags.insert(Flag::DisplacedFromPrior);
                }
            }

            if has_prior && !is_novel && !input.prior_taxa.contains(taxon) {
                row.flags.insert(Flag::NewInDecoration);
            }
        }

        if row.action.is_some() || !row.flags.is_empty() {
            rows.push(row);
        }
    }

    rows.sort_by(|a, b| a.rank.cmp(&b.rank).then_with(|| a.taxon.cmp(&b.taxon)));
    Ok(rows)
}

/// Format an optional value, an empty string if it is not set.
fn optional<T>(value: Option<T>, format: impl Fn(T) -> String) -> String {
    value.map_or_default(format)
}

fn join<T: fmt::Display>(values: impl IntoIterator<Item = T>) -> String {
    values
        .into_iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(LIST_SEPARATOR)
}

/// The report as a TSV, with a header row.
pub(crate) fn to_tsv(rows: &[ReportRow]) -> String {
    let mut out = HEADER.join("\t");
    out.push('\n');
    for row in rows {
        let before = row.before;
        let fields = [
            row.taxon.to_string(),
            row.rank.prefix().to_string(),
            optional(row.action, |action| action.to_string()),
            join(&row.flags),
            optional(before, |b| b.n_genomes.to_string()),
            row.after.n_genomes.to_string(),
            optional(before.and_then(|b| b.bootstrap), |bs| format!("{bs:.1}")),
            optional(row.after.bootstrap, |bs| format!("{bs:.1}")),
            optional(before, |b| format!("{:.3}", b.red)),
            format!("{:.3}", row.after.red),
            join(&row.prior_taxa),
            optional(row.nested_in.as_ref(), ToString::to_string),
            join(&row.sister_taxa),
        ];
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    out
}

/// Write the report to `path`.
pub(crate) fn write_report(path: &Path, rows: &[ReportRow]) -> CuratomaticResult<()> {
    fs::write(path, to_tsv(rows))
        .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))
}

#[cfg(test)]
mod tests;
