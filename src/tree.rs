use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::path::Path;

use ::log::{info, warn};

use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::metadata::Metadata;
use crate::newick::{NewickTree, NodeId};
use crate::node_label::NodeLabel;
use crate::rank::Rank;
use crate::red::{Red, RedValue};
use crate::taxon::{Taxon, generate_unique_rank};
use crate::taxonomy::Taxonomy;

/// Separates the rest of a node label from its RED value.
const RED_SEPARATOR: &str = "|RED=";

/// A species that has no name yet: the rank prefix on its own.
const SPECIES_PLACEHOLDER: &str = "s__";

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

/// A RED decorated tree: leaf nodes are labelled `<gid>|RED=<red>`, and internal nodes
/// `[<bootstrap>][:][<taxonomy>]|RED=<red>`.
pub(crate) struct Tree {
    tree: NewickTree,
    /// The label of every node, indexed by node. Call [`Tree::reload`] after changing the tree.
    node_labels: Vec<NodeLabel>,
    /// The genome ID of every leaf node that is not a zombie.
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
            tree,
            node_labels,
            gids,
            taxon_labels,
        })
    }

    /// Recalculate everything that is derived from the labels of the nodes.
    pub(crate) fn reload(&mut self) -> CuratomaticResult<()> {
        self.node_labels = Self::get_node_labels(&self.tree)?;
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
            .filter(|gid| !gid.is_zombie())
            .cloned()
            .collect()
    }

    /// Every taxon on the tree, which must be on one node only.
    fn get_taxon_labels(node_labels: &[NodeLabel]) -> CuratomaticResult<BTreeSet<Taxon>> {
        let mut out = BTreeSet::new();
        for taxonomy in node_labels.iter().filter_map(|label| label.taxonomy.as_ref()) {
            for taxon in taxonomy.iter_ranks() {
                if !out.insert(taxon.clone()) {
                    return Err(CuratomaticError::NonUniqueTaxon(taxon.to_string()));
                }
            }
        }
        Ok(out)
    }

    fn label(&self, node: NodeId) -> CuratomaticResult<&NodeLabel> {
        self.node_labels
            .get(node.index())
            .ok_or_else(|| CuratomaticError::Invariant(format!("node {node} has no label")))
    }

    /// Read a tree from a Newick file.
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        Self::from_newick(NewickTree::from_file(path)?)
    }

    /// Build the tree from a parsed Newick tree, see [`Tree::from_file`].
    fn from_newick(mut tree: NewickTree) -> CuratomaticResult<Self> {
        // Convert any placeholders (e.g. s__) that exist on the basal node to an unresolved name.
        // Only the placeholder is replaced: any other taxa on the node (e.g. its genus) are kept.
        // TODO: This is a hack, fix it.
        let leaves: Vec<NodeId> = tree.leaves().collect();
        for leaf in leaves {
            let gid = GenomeId::new(
                leaf_label(&tree, leaf)?
                    .split("|RED")
                    .next()
                    .unwrap_or_default(),
            );
            if gid.is_zombie() {
                continue;
            }

            let parent = tree.parent(leaf).ok_or_else(|| {
                CuratomaticError::Invariant(format!("leaf {gid} has no parent node"))
            })?;
            let parent_label = tree.label(parent).unwrap_or_default();
            let (red, bootstrap, taxon) = parse_node_label(parent_label)?;
            if red != RedValue::LEAF || bootstrap.is_some() {
                return Err(CuratomaticError::InvalidNodeLabel(parent_label.to_owned()));
            }

            let taxa: Vec<&str> = taxon.unwrap_or_default().split(';').map(str::trim).collect();
            if taxa.contains(&SPECIES_PLACEHOLDER) {
                let unresolved = format!("{SPECIES_PLACEHOLDER}{{unresolved}} {gid}");
                let taxa = taxa
                    .into_iter()
                    .map(|taxon| if taxon == SPECIES_PLACEHOLDER { unresolved.as_str() } else { taxon })
                    .collect::<Vec<_>>()
                    .join("; ");
                tree.set_label(parent, format!("{taxa}{RED_SEPARATOR}{red:.3}"));
            }
        }

        Self::new(tree)
    }

    /// Write the tree to a Newick file.
    pub(crate) fn to_file(&self, path: &Path) -> CuratomaticResult<()> {
        self.tree.to_file(path)
    }

    /// The number of leaf nodes that are not zombies, and the number that are.
    pub(crate) fn n_nodes_and_zombies(&self) -> CuratomaticResult<(usize, usize)> {
        let (mut n_nodes, mut n_zombies) = (0_usize, 0_usize);
        for leaf in self.tree.leaves() {
            if leaf_gid(&self.tree, leaf).is_zombie() {
                n_zombies = n_zombies.saturating_add(1);
            } else {
                n_nodes = n_nodes.saturating_add(1);
            }
        }
        Ok((n_nodes, n_zombies))
    }

    /// The leaf nodes that are missing one or more ranks, and the ranks that they are missing.
    ///
    /// The ranks are those not assigned on the path from the leaf to the root. Zombies are skipped.
    pub(crate) fn get_leaf_nodes_missing_rank(&self) -> CuratomaticResult<MissingRanks> {
        let mut out = MissingRanks::new();

        for leaf in self.tree.leaves() {
            if leaf_gid(&self.tree, leaf).is_zombie() {
                continue;
            }

            // Go up the tree and read the taxonomy from it.
            let mut taxonomy = Taxonomy::default();
            let mut current = self.tree.parent(leaf);
            while let Some(node) = current {
                if let Some(node_taxonomy) = NodeLabel::from_node(&self.tree, node)?.taxonomy {
                    for taxon in node_taxonomy.iter_ranks() {
                        taxonomy.add_taxon(taxon.clone())?;
                    }
                }

                // Stop early if there is nothing left to find.
                current = if taxonomy.all_assigned() {
                    None
                } else {
                    self.tree.parent(node)
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
        for node in self.tree.ids() {
            let label = self
                .tree
                .label(node)
                .ok_or_else(|| CuratomaticError::InvalidNodeLabel(format!("node {node}")))?;
            let label = label.split('|').next().unwrap_or_default().to_owned();
            self.tree.set_label(node, label);
        }
        Ok(())
    }

    /// Every taxon that is on the tree, read from the labels of the nodes themselves.
    pub(crate) fn extract_all_taxa_from_tree(&self) -> CuratomaticResult<BTreeSet<Taxon>> {
        let mut out = BTreeSet::new();
        for node in self.tree.ids() {
            if let Some(taxonomy) = NodeLabel::from_node(&self.tree, node)?.taxonomy {
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

        for (node, node_label) in self.tree.ids().zip(&self.node_labels) {
            // Skip those that are not set (i.e. leaf nodes).
            let (Some(taxonomy), Some(original_bs)) = (&node_label.taxonomy, node_label.bootstrap)
            else {
                continue;
            };
            let original_red_delta = median_red.distance(node_label.red);

            for taxon in taxonomy.get_rank(rank) {
                // Check up the tree (this is valid for all cases).
                let mut current = self.tree.parent(node);
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

                    current = self.tree.parent(cur_node);
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
        self.tree.descendants(node).any(|desc| {
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
        self.tree.descendants(node).any(|desc| {
            let Some(taxonomy) = self.node_labels.get(desc.index()).and_then(|l| l.taxonomy.as_ref())
            else {
                return false;
            };
            let named = || {
                taxonomy
                    .get_rank(rank)
                    .iter()
                    .filter(|other| other.as_str().chars().count() > 3)
            };
            named().next().is_some() && !named().any(|other| other == taxon)
        })
    }

    /// The taxa of `rank` that the metadata assigns to the genomes below each child of `node`.
    fn child_rank_taxa<'a>(
        &self,
        node: NodeId,
        rank: Rank,
        gid_to_taxonomy: &BTreeMap<&GenomeId, &'a Taxonomy>,
    ) -> CuratomaticResult<BTreeMap<NodeId, BTreeSet<&'a Taxon>>> {
        let mut out: BTreeMap<NodeId, BTreeSet<&'a Taxon>> = BTreeMap::new();
        for &child in self.tree.children(node) {
            for desc in self.tree.descendants(child) {
                if !self.tree.is_leaf(desc) {
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
    /// Call [`Tree::reload`] afterwards.
    pub(crate) fn process_refinement(
        &mut self,
        refinement: &Refinement,
        rank: Rank,
        min_bs: BootstrapValue,
        red_dict: &Red,
    ) -> CuratomaticResult<()> {
        if refinement.entries.is_empty() {
            return Ok(());
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
                    (from_node, (original_red_delta, original_bs, original_bs >= min_bs)),
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
                self.tree.label(original_node).unwrap_or("None"),
                self.tree.label(best_node).unwrap_or("None")
            );
            transfer_taxon(&mut self.tree, original_node, best_node, taxon)?;
        }
        Ok(())
    }
}

/// Parse an internal node label into its RED, bootstrap and taxonomy (if they are set).
pub(crate) fn parse_node_label(
    label: &str,
) -> CuratomaticResult<(RedValue, Option<BootstrapValue>, Option<&str>)> {
    let invalid = || CuratomaticError::InvalidNodeLabel(label.to_owned());
    let float = |value: &str| {
        value
            .parse::<BootstrapValue>()
            .map_err(|_not_a_float| invalid())
    };

    let mut parts = label.trim().split(RED_SEPARATOR);
    let (Some(prefix), Some(suffix), None) = (parts.next(), parts.next(), parts.next()) else {
        return Err(invalid());
    };
    let red = suffix.parse::<RedValue>().map_err(|_not_a_float| invalid())?;

    // Either the taxonomy, the bootstrap, or the bootstrap and taxonomy.
    if prefix.contains(':') {
        let mut parts = prefix.split(':');
        let (Some(bootstrap), Some(taxon), None) = (parts.next(), parts.next(), parts.next())
        else {
            return Err(invalid());
        };
        Ok((red, Some(float(bootstrap)?), Some(taxon)))
    } else if let Ok(bootstrap) = float(prefix) {
        Ok((red, Some(bootstrap), None))
    } else {
        Ok((red, None, Some(prefix)))
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

    // A leaf is labelled with its genome, it cannot hold a taxon.
    if let Some(gid) = &new_info.gid {
        return Err(CuratomaticError::Invariant(format!(
            "{taxon} cannot be moved to the leaf node {gid}"
        )));
    }

    let old_taxonomy = old_info.taxonomy.as_mut().ok_or_else(|| {
        CuratomaticError::Invariant(format!("node {old_node} has no taxonomy to remove {taxon} from"))
    })?;
    old_taxonomy.remove_taxon(taxon)?;

    // A node without a bootstrap is labelled by its taxonomy alone, so removing its last taxon
    // would leave nothing to label it with. Neither label has been changed at this point.
    if old_info.bootstrap.is_none() && old_taxonomy.iter_ranks().next().is_none() {
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
/// Each is placed on the highest node that holds no other taxon of the rank below it, has enough
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
            let mut current = tree.tree.parent(leaf);
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
                current = tree.tree.parent(node);
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
            let &candidate = candidates.last().ok_or_else(|| CuratomaticError::NoCandidateNode {
                rank: rank.prefix(),
                gid: leaf_gid(&tree.tree, *leaf),
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
                .map(|&leaf| NodeLabel::from_node(&tree.tree, leaf).map(|label| label.gid))
                .collect::<CuratomaticResult<Vec<_>>>()?
                .into_iter()
                .flatten()
                .collect::<BTreeSet<GenomeId>>();
            let new_taxon = generate_unique_rank(*rank, &gids, metadata, &known_taxa)?;
            known_taxa.insert(new_taxon.clone());
            new_taxa_created.insert(new_taxon.clone());

            // Add the taxon to the label of the node.
            let mut info = NodeLabel::from_node(&tree.tree, *candidate)?;
            info.taxonomy
                .get_or_insert_with(Taxonomy::default)
                .add_taxon(new_taxon)?;
            tree.tree.set_label(*candidate, info.to_label()?);
        }
    }

    Ok(new_taxa_created)
}

/// The label of a leaf node.
fn leaf_label(tree: &NewickTree, leaf: NodeId) -> CuratomaticResult<&str> {
    tree.label(leaf)
        .ok_or_else(|| CuratomaticError::InvalidNodeLabel(format!("leaf {leaf} has no label")))
}

/// The genome ID of a leaf node, which is that of a zombie if the leaf is one.
fn leaf_gid(tree: &NewickTree, leaf: NodeId) -> GenomeId {
    GenomeId::new(
        tree.label(leaf)
            .and_then(|label| label.split(RED_SEPARATOR).next())
            .unwrap_or_default(),
    )
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
