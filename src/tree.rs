use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::path::Path;

use ::log::{info, warn};
use strum::IntoEnumIterator;

use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::metadata::Metadata;
use crate::newick::{NewickTree, NodeId};
use crate::node_label::{NodeLabel, strip_red};
use crate::rank::Rank;
use crate::red::Red;
use crate::taxon::{Taxon, generate_unique_rank};
use crate::taxonomy::Taxonomy;
use crate::taxonomy_file::TaxonomyFile;

/// The ranks that each leaf node is missing, in the order that the leaves occur in the tree.
pub(crate) type MissingRanks = BTreeMap<NodeId, BTreeSet<Rank>>;

/// Where each taxon could be moved to, in the order they were found.
///
/// Each entry is a taxon, and its candidate moves: `(to_node, from_node)`.
#[derive(Debug, Default)]
pub(crate) struct Refinement {
    entries: Vec<(Taxon, Vec<(NodeId, NodeId)>)>,
    index: HashMap<Taxon, usize>,
}

impl Refinement {
    fn push(&mut self, taxon: &Taxon, to_node: NodeId, from_node: NodeId) {
        let idx = *self.index.entry(taxon.clone()).or_insert_with(|| {
            self.entries.push((taxon.clone(), Vec::new()));
            self.entries.len().saturating_sub(1)
        });
        if let Some((_, moves)) = self.entries.get_mut(idx) {
            moves.push((to_node, from_node));
        }
    }
}

/// What [`Tree::assign_leaf_taxonomy`] did.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct LeafAssignment {
    /// The number of leaves that were given one or more taxa.
    pub(crate) n_leaves: usize,
    /// The number of taxa that were assigned to a leaf, for each rank.
    pub(crate) n_taxa: BTreeMap<Rank, usize>,
    /// Taxa that were not assigned: they are already on the tree, or more than one leaf would get
    /// them. Each taxon is on one node only.
    pub(crate) skipped: BTreeSet<Taxon>,
}

/// A taxon that [`Tree::process_refinement`] moved from one node to another.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Move {
    pub(crate) taxon: Taxon,
    pub(crate) from: NodeId,
    pub(crate) to: NodeId,
}

/// A RED decorated tree: leaf nodes are labelled `<gid>|RED=<red>`, and internal nodes
/// `[<bootstrap>][:][<taxonomy>]|RED=<red>`.
pub(crate) struct Tree {
    newick: NewickTree,
    /// The label of every node, indexed by node. Call [`Tree::reload`] after changing the tree.
    node_labels: Vec<NodeLabel>,
    /// The genome ID of every leaf node.
    pub(crate) gids: BTreeSet<GenomeId>,
    /// Every taxon that is on the tree.
    pub(crate) taxon_labels: BTreeSet<Taxon>,
}

impl Tree {
    fn new(tree: NewickTree) -> CuratomaticResult<Self> {
        let node_labels = Self::get_node_labels(&tree)?;
        let gids = Self::get_gids(&tree, &node_labels);
        let taxon_labels = Self::get_taxon_labels(&node_labels)?;
        Ok(Tree {
            newick: tree,
            node_labels,
            gids,
            taxon_labels,
        })
    }

    /// Recalculate everything that is derived from the labels of the nodes.
    pub(crate) fn reload(&mut self) -> CuratomaticResult<()> {
        self.node_labels = Self::get_node_labels(&self.newick)?;
        self.taxon_labels = Self::get_taxon_labels(&self.node_labels)?;
        Ok(())
    }

    fn get_node_labels(tree: &NewickTree) -> CuratomaticResult<Vec<NodeLabel>> {
        tree.ids()
            .map(|node| NodeLabel::from_node(tree, node))
            .collect()
    }

    fn get_gids(tree: &NewickTree, node_labels: &[NodeLabel]) -> BTreeSet<GenomeId> {
        tree.leaves()
            .filter_map(|node| node_labels.get(node.index())?.gid.as_ref())
            .cloned()
            .collect()
    }

    /// Every taxon on the tree, which must be on one node only (except [`Taxon::unresolvable`]).
    fn get_taxon_labels(node_labels: &[NodeLabel]) -> CuratomaticResult<BTreeSet<Taxon>> {
        let mut out = BTreeSet::new();
        for taxonomy in node_labels
            .iter()
            .filter_map(|label| label.taxonomy.as_ref())
        {
            for taxon in taxonomy.iter_ranks() {
                if !out.insert(taxon.clone()) && !taxon.is_unresolvable() {
                    return Err(CuratomaticError::NonUniqueTaxon(taxon.to_string()));
                }
            }
        }
        Ok(out)
    }

    pub(crate) fn label(&self, node: NodeId) -> CuratomaticResult<&NodeLabel> {
        self.node_labels
            .get(node.index())
            .ok_or_else(|| CuratomaticError::Invariant(format!("node {node} has no label")))
    }

    /// The underlying Newick tree.
    pub(crate) fn newick(&self) -> &NewickTree {
        &self.newick
    }

    /// The node that holds each taxon. Call [`Tree::reload`] first if the tree was changed.
    pub(crate) fn taxon_nodes(&self) -> BTreeMap<Taxon, NodeId> {
        self.newick
            .ids()
            .zip(&self.node_labels)
            .filter_map(|(node, label)| Some((node, label.taxonomy.as_ref()?)))
            .flat_map(|(node, taxonomy)| {
                taxonomy
                    .iter_ranks()
                    .map(move |taxon| (taxon.clone(), node))
            })
            .collect()
    }

    /// Read a tree from a Newick file.
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        Self::from_newick(NewickTree::from_file(path)?)
    }

    /// Build the tree from a parsed Newick tree, see [`Tree::from_file`].
    pub(crate) fn from_newick(tree: NewickTree) -> CuratomaticResult<Self> {
        Self::new(tree)
    }

    /// Write the tree to a Newick file.
    pub(crate) fn to_file(&self, path: &Path) -> CuratomaticResult<()> {
        self.newick.to_file(path)
    }

    /// The number of leaf nodes.
    pub(crate) fn n_leaves(&self) -> usize {
        self.newick.leaves().count()
    }

    /// The genome ID of a leaf, `None` if `node` is not a leaf.
    fn leaf_gid(&self, node: NodeId) -> Option<&GenomeId> {
        self.node_labels.get(node.index())?.gid.as_ref()
    }

    /// The lowest rank of any taxon (placeholders included) on the path from each node (included)
    /// to the root, indexed by node.
    fn lowest_rank_on_path(&self) -> CuratomaticResult<Vec<Option<Rank>>> {
        let mut out: Vec<Option<Rank>> = vec![None; self.node_labels.len()];
        // Post-order puts every parent after its children, so in reverse the parent comes first.
        let nodes: Vec<NodeId> = self.newick.ids().collect();
        for &node in nodes.iter().rev() {
            let mut lowest = self
                .newick
                .parent(node)
                .and_then(|parent| out.get(parent.index()).copied().flatten());
            if let Some(taxonomy) = &self.label(node)?.taxonomy {
                for taxon in taxonomy.iter_ranks() {
                    lowest = lowest.max(Some(taxon.rank()?));
                }
            }
            if let Some(slot) = out.get_mut(node.index()) {
                *slot = lowest;
            }
        }
        Ok(out)
    }

    /// Give each leaf the ranks that the tree can't hold, from the taxonomy file.
    ///
    /// Each leaf gets every named taxon of the file whose rank is below the
    /// lowest rank on its path to the root, e.g. only its species if the lowest is a genus. A gap
    /// above the lowest rank is never filled, and placeholders (e.g. `g__`) are never assigned.
    /// The taxa already on the tree are never changed. A taxon that is already on the tree, or
    /// that more than one leaf would get, is not assigned (see [`LeafAssignment::skipped`]).
    ///
    /// A leaf whose genome is not in the file is an error. Reloads the tree.
    pub(crate) fn assign_leaf_taxonomy(
        &mut self,
        file: &TaxonomyFile,
    ) -> CuratomaticResult<LeafAssignment> {
        let lowest = self.lowest_rank_on_path()?;

        // The taxa that each leaf would get, and how many leaves would get each taxon.
        let mut wanted: Vec<(NodeId, Vec<&Taxon>)> = Vec::new();
        let mut n_wanting: HashMap<&Taxon, usize> = HashMap::new();
        for leaf in self.newick.leaves() {
            let Some(gid) = self.leaf_gid(leaf) else {
                continue;
            };
            let taxonomy =
                file.get(gid)
                    .ok_or_else(|| CuratomaticError::MissingFromTaxonomyFile {
                        n: 1,
                        first: gid.clone(),
                    })?;
            let below = lowest.get(leaf.index()).copied().flatten();
            let mut taxa = Vec::new();
            for taxon in taxonomy.iter_ranks().filter(|taxon| taxon.is_named()) {
                if below.is_none_or(|below| taxon.rank().is_ok_and(|rank| rank > below)) {
                    taxa.push(taxon);
                    let n = n_wanting.entry(taxon).or_default();
                    *n = n.saturating_add(1);
                }
            }
            if !taxa.is_empty() {
                wanted.push((leaf, taxa));
            }
        }

        let mut out = LeafAssignment::default();
        for (leaf, taxa) in wanted {
            let mut info = self.label(leaf)?.clone();
            let mut assigned = false;
            for taxon in taxa {
                if self.taxon_labels.contains(taxon) || n_wanting.get(taxon).is_some_and(|&n| n > 1)
                {
                    out.skipped.insert(taxon.clone());
                    continue;
                }
                let n = out.n_taxa.entry(taxon.rank()?).or_default();
                *n = n.saturating_add(1);
                info.taxonomy
                    .get_or_insert_with(Taxonomy::default)
                    .add_taxon(taxon.clone())?;
                assigned = true;
            }
            if assigned {
                out.n_leaves = out.n_leaves.saturating_add(1);
                self.newick.set_label(leaf, info.to_label()?);
            }
        }

        self.reload()?;
        Ok(out)
    }

    /// The taxonomy of each leaf, in tree order.
    ///
    /// The taxonomy is read from the taxa on the leaf and every ancestor, one taxon for each rank
    /// from the domain to the species. When a rank is on more than one node of the path (nested
    /// taxa) the one nearest the leaf is used, and when one node holds more than one, the first.
    /// Placeholders are ignored, a rank without a named taxon is given as its placeholder (e.g.
    /// `g__`). Call [`Tree::reload`] first if the tree was changed.
    pub(crate) fn leaf_taxonomies(&self) -> CuratomaticResult<Vec<(GenomeId, Vec<Taxon>)>> {
        let mut out = Vec::new();
        for leaf in self.newick.leaves() {
            let Some(gid) = self.leaf_gid(leaf) else {
                continue;
            };
            let mut named: BTreeMap<Rank, &Taxon> = BTreeMap::new();
            let mut current = Some(leaf);
            while let Some(node) = current {
                if let Some(taxonomy) = &self.label(node)?.taxonomy {
                    for taxon in taxonomy.iter_ranks().filter(|taxon| taxon.is_named()) {
                        named.entry(taxon.rank()?).or_insert(taxon);
                    }
                }
                current = self.newick.parent(node);
            }
            let taxa = Rank::iter()
                .map(|rank| {
                    named.get(&rank).map_or_else(
                        || Taxon::new(format!("{}__", rank.prefix())),
                        |&taxon| taxon.clone(),
                    )
                })
                .collect();
            out.push((gid.clone(), taxa));
        }
        Ok(out)
    }

    /// The leaf nodes that are missing one or more ranks, and the ranks that they are missing.
    ///
    /// The ranks are those not assigned on the path from the leaf (included) to the root. Call [`Tree::reload`] first if the tree was changed.
    pub(crate) fn get_leaf_nodes_missing_rank(&self) -> CuratomaticResult<MissingRanks> {
        let mut out = MissingRanks::new();

        for leaf in self.newick.leaves() {
            // Go up the tree and read the taxonomy from it.
            let mut taxonomy = Taxonomy::default();
            let mut current = Some(leaf);
            while let Some(node) = current {
                if let Some(node_taxonomy) = &self.label(node)?.taxonomy {
                    for taxon in node_taxonomy.iter_ranks() {
                        taxonomy.add_taxon(taxon.clone())?;
                    }
                }

                // Stop early if there is nothing left to find.
                current = if taxonomy.all_assigned() {
                    None
                } else {
                    self.newick.parent(node)
                };
            }

            if !taxonomy.all_assigned() {
                out.insert(leaf, taxonomy.missing_ranks());
            }
        }
        Ok(out)
    }

    /// Remove the RED value from the label of every node.
    pub(crate) fn remove_red_labels(&mut self) -> CuratomaticResult<()> {
        for node in self.newick.ids() {
            let label = self
                .newick
                .label(node)
                .ok_or_else(|| CuratomaticError::InvalidNodeLabel(format!("node {node}")))?;
            let label = strip_red(label)?.to_owned();
            self.newick.set_label(node, label);
        }
        Ok(())
    }

    /// Every taxon that is on the tree, read from the labels of the nodes themselves.
    pub(crate) fn extract_all_taxa_from_tree(&self) -> CuratomaticResult<BTreeSet<Taxon>> {
        let mut out = BTreeSet::new();
        for node in self.newick.ids() {
            if let Some(taxonomy) = NodeLabel::from_node(&self.newick, node)?.taxonomy {
                out.extend(taxonomy.iter_ranks().cloned());
            }
        }
        Ok(out)
    }

    /// Find the nodes that each taxon of `rank` could be moved to.
    ///
    /// A taxon can be moved up the tree while its RED stays within the range of the rank, and
    /// down the tree towards the genomes that the metadata assigns it to. Designed to be called
    /// from phylum to genus.
    pub(crate) fn refine_existing_labels(
        &self,
        rank: Rank,
        min_bs: BootstrapValue,
        red_dict: &Red,
        metadata: &Metadata,
    ) -> CuratomaticResult<Refinement> {
        let mut out = Refinement::default();
        let gid_to_taxonomy = metadata.get_taxonomy_from_gids(&self.gids);

        let child_rank = rank.child().ok_or_else(|| {
            CuratomaticError::Invariant("the species cannot be refined, it has no child".to_owned())
        })?;
        let median_red = red_dict.get(rank)?;
        let parent_red_thresh = red_dict.get(parent_rank(rank)?)?.midpoint(median_red);
        let child_red_thresh = red_dict.get(child_rank)?.midpoint(median_red);

        for (node, node_label) in self.newick.ids().zip(&self.node_labels) {
            // Skip nodes without a bootstrap (i.e. leaf nodes), which are never refined.
            let (Some(taxonomy), Some(original_bs)) = (&node_label.taxonomy, node_label.bootstrap)
            else {
                continue;
            };
            let original_red_delta = median_red.distance(node_label.red);

            for taxon in taxonomy.get_rank(rank) {
                // Check up the tree (this is valid for all cases). A taxon nested inside another
                // taxon of the same rank is not moved up: it would climb onto or past the node of
                // that taxon, taking its genomes.
                let mut current = if self.nested_inside(node, rank, taxon).is_some() {
                    None
                } else {
                    self.newick.parent(node)
                };
                while let Some(cur_node) = current {
                    let cur_label = self.label(cur_node)?;
                    let cur_red = cur_label.red;

                    // Stop once the maximum RED range is exceeded.
                    if cur_red < parent_red_thresh {
                        break;
                    }

                    // Check if this node is within the RED range.
                    if cur_red <= child_red_thresh {
                        // Make sure that this will actually improve the bootstrap support.
                        //
                        // Every internal node is expected to have a bootstrap value. The only
                        // nodes without one are the species nodes directly above the leaves,
                        // which are never an ancestor of a node that holds a bootstrap. An
                        // ancestor without one (e.g. a root labelled `d__Bacteria|RED=0.000`)
                        // is a malformed tree, and is reported as an error rather than guessed at.
                        let cur_bs = required_bootstrap(cur_label)?;
                        if cur_bs >= min_bs || cur_bs > original_bs {
                            // Stop if this would incorporate any nested ranks.
                            if self.would_cause_conflicts(cur_node, rank, taxon) {
                                break;
                            }

                            if median_red.distance(cur_red) < original_red_delta {
                                // This brings the RED closer to the median.
                                out.push(taxon, cur_node, node);
                            } else if original_bs < min_bs {
                                // This moves it further from the median, only do so to
                                // improve on a node that is not supported.
                                warn!("UNTESTED CODE");
                                out.push(taxon, cur_node, node);
                            }
                        }
                    }

                    current = self.newick.parent(cur_node);
                }

                // A well supported taxon is not moved down if it would move the RED further
                // away from the median.
                if original_bs >= min_bs && node_label.red > child_red_thresh {
                    continue;
                }

                // Find which children contain genomes that the metadata assigns to the taxon.
                let child_taxa = self.child_rank_taxa(node, rank, &gid_to_taxonomy)?;

                // Moving the label down would cause polyphyly, which is not handled.
                if polyphyly(&child_taxa, taxon) {
                    continue;
                }

                // Otherwise, go over the children that wouldn't cause polyphyly.
                let mut queue: VecDeque<NodeId> = containing(&child_taxa, taxon).collect();
                while let Some(cur_node) = queue.pop_front() {
                    let cur_label = self.label(cur_node)?;

                    match cur_label.bootstrap {
                        // A leaf node is included regardless of the RED range, given that the
                        // prior node had bad support. Do not go down any further.
                        None => {
                            out.push(taxon, cur_node, node);
                            continue;
                        }
                        Some(cur_bs) if cur_bs >= min_bs => {
                            out.push(taxon, cur_node, node);

                            // Given that a good node was found, stop if straying further from
                            // the median.
                            if median_red.distance(cur_label.red) > original_red_delta {
                                continue;
                            }
                        }
                        // This is an improvement, keep searching.
                        Some(cur_bs) if cur_bs > original_bs => {
                            out.push(taxon, cur_node, node);
                        }
                        Some(_) => {}
                    }

                    let child_taxa = self.child_rank_taxa(cur_node, rank, &gid_to_taxonomy)?;
                    if polyphyly(&child_taxa, taxon) {
                        continue;
                    }
                    queue.extend(containing(&child_taxa, taxon));
                }
            }
        }

        Ok(out)
    }

    /// True if any node below `node` has a taxon for `rank`.
    fn any_descendant_has_rank(&self, node: NodeId, rank: Rank) -> bool {
        self.newick.descendants(node).any(|desc| {
            self.node_labels
                .get(desc.index())
                .and_then(|label| label.taxonomy.as_ref())
                .is_some_and(|taxonomy| !taxonomy.get_rank(rank).is_empty())
        })
    }

    /// True if a node below `node` has a different taxon for `rank` than `taxon`.
    ///
    /// Taxa that are only the rank prefix (e.g. `g__`) are placeholders, and are ignored.
    fn would_cause_conflicts(&self, node: NodeId, rank: Rank, taxon: &Taxon) -> bool {
        self.newick.descendants(node).any(|desc| {
            let Some(taxonomy) = self
                .node_labels
                .get(desc.index())
                .and_then(|l| l.taxonomy.as_ref())
            else {
                return false;
            };
            let named = || {
                taxonomy
                    .get_rank(rank)
                    .iter()
                    .filter(|other| other.is_named())
            };
            named().next().is_some() && !named().any(|other| other == taxon)
        })
    }

    /// The closest named taxon of `rank`, other than `taxon`, on a node above `node`.
    ///
    /// If there is one, then `taxon` is nested inside it. Placeholders (e.g. `g__`) are ignored.
    pub(crate) fn nested_inside(&self, node: NodeId, rank: Rank, taxon: &Taxon) -> Option<&Taxon> {
        let mut current = self.newick.parent(node);
        while let Some(ancestor) = current {
            let found = self
                .node_labels
                .get(ancestor.index())
                .and_then(|label| label.taxonomy.as_ref())
                .and_then(|taxonomy| {
                    taxonomy
                        .get_rank(rank)
                        .iter()
                        .find(|other| other.is_named() && *other != taxon)
                });
            if found.is_some() {
                return found;
            }
            current = self.newick.parent(ancestor);
        }
        None
    }

    /// The taxa of `rank` that the metadata assigns to the genomes at or below each child of
    /// `node` (a child may be a leaf itself).
    fn child_rank_taxa<'a>(
        &self,
        node: NodeId,
        rank: Rank,
        gid_to_taxonomy: &BTreeMap<&GenomeId, &'a Taxonomy>,
    ) -> CuratomaticResult<BTreeMap<NodeId, BTreeSet<&'a Taxon>>> {
        let mut out: BTreeMap<NodeId, BTreeSet<&'a Taxon>> = BTreeMap::new();
        for &child in self.newick.children(node) {
            for desc in self.newick.descendants(child).chain([child]) {
                if !self.newick.is_leaf(desc) {
                    continue;
                }
                let gid = self.label(desc)?.gid.as_ref();
                if let Some(taxonomy) = gid.and_then(|gid| gid_to_taxonomy.get(gid)) {
                    out.entry(child)
                        .or_default()
                        .extend(taxonomy.get_rank(rank));
                }
            }
        }
        Ok(out)
    }

    /// Apply the best of the candidate moves found by [`Tree::refine_existing_labels`].
    ///
    /// The candidates for a taxon are sorted by:
    /// 1. Whether the bootstrap support meets the minimum.
    /// 2. The RED delta to the median of the rank.
    /// 3. The bootstrap support, descending.
    ///
    /// Returns the moves that were made. Call [`Tree::reload`] afterwards.
    pub(crate) fn process_refinement(
        &mut self,
        refinement: &Refinement,
        rank: Rank,
        min_bs: BootstrapValue,
        red_dict: &Red,
    ) -> CuratomaticResult<Vec<Move>> {
        let mut moves = Vec::new();
        if refinement.entries.is_empty() {
            return Ok(moves);
        }
        let median_red = red_dict.get(rank)?;

        for (taxon, node_changes) in &refinement.entries {
            let Some(&(_, original_node)) = node_changes.first() else {
                continue;
            };
            let original_red_delta = median_red.distance(self.label(original_node)?.red);

            // The (RED delta, bootstrap, supported) of each node, in the order they were added.
            let mut changes: Vec<(NodeId, (f64, BootstrapValue, bool))> = Vec::new();
            for &(to_node, from_node) in node_changes {
                let target_label = self.label(to_node)?;
                let target_red_delta = median_red.distance(target_label.red);
                let target_bs = target_label.bootstrap.unwrap_or(BootstrapValue::LEAF);
                let original_bs = required_bootstrap(self.label(from_node)?)?;

                for (node, change) in [
                    (to_node, (target_red_delta, target_bs, target_bs >= min_bs)),
                    (
                        from_node,
                        (original_red_delta, original_bs, original_bs >= min_bs),
                    ),
                ] {
                    match changes.iter_mut().find(|(existing, _)| *existing == node) {
                        Some((_, existing)) => *existing = change,
                        None => changes.push((node, change)),
                    }
                }
            }

            changes.sort_by(|(_, a), (_, b)| {
                b.2.cmp(&a.2)
                    .then(a.0.total_cmp(&b.0))
                    .then(b.1.total_cmp(a.1))
            });
            let Some(&(best_node, _)) = changes.first() else {
                continue;
            };

            // Nothing to do here, ideally this wouldn't happen.
            if best_node == original_node {
                continue;
            }

            info!(
                "Moving {taxon} from {} to {}",
                self.newick.label(original_node).unwrap_or("None"),
                self.newick.label(best_node).unwrap_or("None")
            );
            transfer_taxon(&mut self.newick, original_node, best_node, taxon)?;
            moves.push(Move {
                taxon: taxon.clone(),
                from: original_node,
                to: best_node,
            });
        }
        Ok(moves)
    }
}

/// Move a taxon from one node to another, rewriting the labels of both.
fn transfer_taxon(
    tree: &mut NewickTree,
    old_node: NodeId,
    new_node: NodeId,
    taxon: &Taxon,
) -> CuratomaticResult<()> {
    let mut old_info = NodeLabel::from_node(tree, old_node)?;
    let mut new_info = NodeLabel::from_node(tree, new_node)?;

    let old_taxonomy = old_info.taxonomy.as_mut().ok_or_else(|| {
        CuratomaticError::Invariant(format!(
            "node {old_node} has no taxonomy to remove {taxon} from"
        ))
    })?;
    old_taxonomy.remove_taxon(taxon)?;

    // An internal node without a bootstrap is labelled by its taxonomy alone, so removing its last
    // taxon would leave nothing to label it with (a leaf keeps its genome ID). Neither label has
    // been changed at this point.
    if old_info.gid.is_none()
        && old_info.bootstrap.is_none()
        && old_taxonomy.iter_ranks().next().is_none()
    {
        return Err(CuratomaticError::Invariant(format!(
            "moving {taxon} would leave node {old_node} with neither a bootstrap nor a taxonomy"
        )));
    }
    new_info
        .taxonomy
        .get_or_insert_with(Taxonomy::default)
        .add_taxon(taxon.clone())?;

    tree.set_label(old_node, old_info.to_label()?);
    tree.set_label(new_node, new_info.to_label()?);
    Ok(())
}

/// Create a novel taxon for every rank that a leaf is missing, and place them on the tree.
///
/// Each is placed on the highest node (the leaf included) that holds no other taxon of the rank
/// below it, has enough
/// bootstrap support, and has a RED that is past the midpoint between the rank and its parent.
/// Leaves that share the highest node share the taxon. The ranks are filled from highest to
/// lowest. `known_taxa` are the taxa that already exist, and are not reused. Returns the taxa
/// that were created.
///
/// The tree must be up to date: call [`Tree::reload`] first.
pub(crate) fn fill_in_missing_ranks(
    tree: &mut Tree,
    leaf_to_ranks_missing: &MissingRanks,
    min_bs: BootstrapValue,
    red_dict: &Red,
    metadata: &Metadata,
    mut known_taxa: BTreeSet<Taxon>,
) -> CuratomaticResult<BTreeSet<Taxon>> {
    // For each rank, the candidate nodes (lowest to highest) of each leaf missing it.
    let mut rank_to_leaf_candidates: BTreeMap<Rank, Vec<(NodeId, Vec<NodeId>)>> = BTreeMap::new();

    for (&leaf, missing_ranks) in leaf_to_ranks_missing {
        // Try find the highest node placement for each rank descending.
        for &rank in missing_ranks {
            // A novel domain cannot be placed, as there is no parent rank to bound its RED.
            let parent_red = red_dict.get(parent_rank(rank)?)?;
            let midpoint_red = parent_red.halfway_to(red_dict.get(rank)?);

            // Go up the tree to find a suitable node that won't conflict.
            let mut candidates = Vec::new();
            let mut current = Some(leaf);
            while let Some(node) = current {
                if tree.any_descendant_has_rank(node, rank) {
                    break;
                }

                // Check RED and bootstrap support.
                let label = tree.label(node)?;
                let bs_supported = label.bootstrap.is_none_or(|bs| bs >= min_bs);
                if bs_supported && label.red >= midpoint_red {
                    candidates.push(node);
                }
                current = tree.newick.parent(node);
            }

            rank_to_leaf_candidates
                .entry(rank)
                .or_default()
                .push((leaf, candidates));
        }
    }

    // Process higher to lower ranks.
    let mut new_taxa_created = BTreeSet::new();
    for (rank, leaf_candidates) in &rank_to_leaf_candidates {
        // Leaves share a taxon if they share the highest candidate node.
        let mut candidate_to_leaves: Vec<(NodeId, BTreeSet<NodeId>)> = Vec::new();
        let mut index: HashMap<NodeId, usize> = HashMap::new();
        for (leaf, candidates) in leaf_candidates {
            let &candidate =
                candidates
                    .last()
                    .ok_or_else(|| CuratomaticError::NoCandidateNode {
                        rank: rank.prefix(),
                        gid: tree
                            .leaf_gid(*leaf)
                            .cloned()
                            .unwrap_or_else(|| GenomeId::new(leaf.to_string())),
                    })?;
            let idx = *index.entry(candidate).or_insert_with(|| {
                candidate_to_leaves.push((candidate, BTreeSet::new()));
                candidate_to_leaves.len().saturating_sub(1)
            });
            if let Some((_, leaves)) = candidate_to_leaves.get_mut(idx) {
                leaves.insert(*leaf);
            }
        }

        // Create a novel group for the highest candidate node that covers the descendant taxa.
        for (candidate, leaves) in &candidate_to_leaves {
            let gids = leaves
                .iter()
                .map(|&leaf| NodeLabel::from_node(&tree.newick, leaf).map(|label| label.gid))
                .collect::<CuratomaticResult<Vec<_>>>()?
                .into_iter()
                .flatten()
                .collect::<BTreeSet<GenomeId>>();
            let new_taxon = generate_unique_rank(*rank, &gids, metadata, &known_taxa)?;
            known_taxa.insert(new_taxon.clone());
            new_taxa_created.insert(new_taxon.clone());

            // Add the taxon to the label of the node.
            let mut info = NodeLabel::from_node(&tree.newick, *candidate)?;
            info.taxonomy
                .get_or_insert_with(Taxonomy::default)
                .add_taxon(new_taxon)?;
            tree.newick.set_label(*candidate, info.to_label()?);
        }
    }

    Ok(new_taxa_created)
}

/// The rank above `rank`, which every rank but the domain has.
fn parent_rank(rank: Rank) -> CuratomaticResult<Rank> {
    rank.parent()
        .ok_or(CuratomaticError::NoParentRank(rank.prefix()))
}

/// The bootstrap support of a node, which an internal node must have to be compared.
fn required_bootstrap(label: &NodeLabel) -> CuratomaticResult<BootstrapValue> {
    label.bootstrap.ok_or_else(|| {
        CuratomaticError::Invariant("expected a node with bootstrap support".to_owned())
    })
}

/// The children that contain `taxon`.
fn containing<'a>(
    child_taxa: &'a BTreeMap<NodeId, BTreeSet<&Taxon>>,
    taxon: &'a Taxon,
) -> impl Iterator<Item = NodeId> + 'a {
    child_taxa
        .iter()
        .filter(move |(_, taxa)| taxa.contains(taxon))
        .map(|(&child, _)| child)
}

/// True if more than one child contains `taxon`.
fn polyphyly(child_taxa: &BTreeMap<NodeId, BTreeSet<&Taxon>>, taxon: &Taxon) -> bool {
    containing(child_taxa, taxon).nth(1).is_some()
}

#[cfg(test)]
mod tests;
