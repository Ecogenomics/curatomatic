use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::path::Path;

use ::log::info;
use strum::IntoEnumIterator;

use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::metadata::Metadata;
use crate::newick::{NewickTree, NodeId, escape_label};
use crate::node_label::{NodeLabel, strip_red};
use crate::rank::Rank;
use crate::red::Red;
use crate::taxon::{Taxon, TaxonCounts, generate_unique_rank};
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
}

/// A taxon that [`Tree::process_refinement`] moved from one node to another.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Move {
    pub(crate) taxon: Taxon,
    pub(crate) from: NodeId,
    pub(crate) to: NodeId,
}

/// A leaf of the tree, which is a genome.
///
/// The label of a leaf on the Newick tree is its genome ID, and is never changed: the taxa that the
/// leaf holds are kept here, and are written to the taxonomy file (see [`Tree::leaf_taxonomies`]).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Leaf {
    pub(crate) gid: GenomeId,
    /// The taxa on the leaf itself (not those of its ancestors).
    pub(crate) taxonomy: Taxonomy,
}

/// A RED decorated tree: leaf nodes are labelled `<gid>[|RED=<red>]`, and internal nodes
/// `[<bootstrap>][:][<taxonomy>]|RED=<red>`.
///
/// The taxa of an internal node are kept in its label, and those of a leaf on its [`Leaf`]. Leaf
/// labels are never changed.
pub(crate) struct Tree {
    newick: NewickTree,
    /// Every leaf node, and the taxa it holds.
    leaves: BTreeMap<NodeId, Leaf>,
    /// The label of every node, indexed by node, where the taxonomy of a leaf is that of its
    /// [`Leaf`]. Call [`Tree::reload`] after changing the tree.
    node_labels: Vec<NodeLabel>,
    /// The genome ID of every leaf node.
    pub(crate) gids: BTreeSet<GenomeId>,
    /// Every taxon that is on the tree.
    pub(crate) taxon_labels: BTreeSet<Taxon>,
}

impl Tree {
    fn new(tree: NewickTree) -> CuratomaticResult<Self> {
        let mut leaves = BTreeMap::new();
        for node in tree.leaves() {
            let gid = NodeLabel::from_node(&tree, node)?.gid.ok_or_else(|| {
                CuratomaticError::Invariant(format!("leaf {node} has no genome ID"))
            })?;
            leaves.insert(
                node,
                Leaf {
                    gid,
                    taxonomy: Taxonomy::default(),
                },
            );
        }
        let node_labels = Self::get_node_labels(&tree, &leaves)?;
        let gids = leaves.values().map(|leaf| leaf.gid.clone()).collect();
        let taxon_labels = Self::get_taxon_labels(&node_labels)?;
        Ok(Tree {
            newick: tree,
            leaves,
            node_labels,
            gids,
            taxon_labels,
        })
    }

    /// Recalculate everything that is derived from the labels of the nodes and the leaves.
    pub(crate) fn reload(&mut self) -> CuratomaticResult<()> {
        self.node_labels = Self::get_node_labels(&self.newick, &self.leaves)?;
        self.taxon_labels = Self::get_taxon_labels(&self.node_labels)?;
        Ok(())
    }

    fn get_node_labels(
        tree: &NewickTree,
        leaves: &BTreeMap<NodeId, Leaf>,
    ) -> CuratomaticResult<Vec<NodeLabel>> {
        tree.ids()
            .map(|node| {
                let mut label = NodeLabel::from_node(tree, node)?;
                if let Some(leaf) = leaves.get(&node) {
                    label.taxonomy = non_empty(leaf.taxonomy.clone());
                }
                Ok(label)
            })
            .collect()
    }

    /// The taxa held by a node, read from the tree as it is now (not from [`Tree::label`], which
    /// is only up to date after [`Tree::reload`]).
    fn node_taxonomy(&self, node: NodeId) -> CuratomaticResult<Taxonomy> {
        if let Some(leaf) = self.leaves.get(&node) {
            return Ok(leaf.taxonomy.clone());
        }
        Ok(NodeLabel::from_node(&self.newick, node)?
            .taxonomy
            .unwrap_or_default())
    }

    /// Replace the taxa held by a node: on its [`Leaf`] for a leaf, and in the label of an
    /// internal node. The label of a leaf is never changed.
    fn set_node_taxonomy(&mut self, node: NodeId, taxonomy: Taxonomy) -> CuratomaticResult<()> {
        if let Some(leaf) = self.leaves.get_mut(&node) {
            leaf.taxonomy = taxonomy;
            return Ok(());
        }
        let mut info = NodeLabel::from_node(&self.newick, node)?;
        info.taxonomy = non_empty(taxonomy);
        self.newick.set_internal_label(node, info.to_label()?)
    }

    /// Add a taxon to a node, see [`Tree::set_node_taxonomy`].
    fn add_node_taxon(&mut self, node: NodeId, taxon: Taxon) -> CuratomaticResult<()> {
        let mut taxonomy = self.node_taxonomy(node)?;
        taxonomy.add_taxon(taxon)?;
        self.set_node_taxonomy(node, taxonomy)
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

    /// Every taxon on the tree and the node that holds it, in node order. A taxon on more than one
    /// node (only [`Taxon::unresolvable`]) is listed once for each. Call [`Tree::reload`] first if
    /// the tree was changed.
    pub(crate) fn taxon_nodes(&self) -> Vec<(Taxon, NodeId)> {
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
        self.leaves.get(&node).map(|leaf| &leaf.gid)
    }

    /// The lowest rank of any valid name (see [`Taxon::is_valid_name`]) on each node or below it,
    /// indexed by node, species excluded. `None` if there is no such name at or below the node.
    ///
    /// Species are excluded because every genome has a named species, so counting them would
    /// stop every novel genus.
    fn lowest_valid_rank_below(&self) -> CuratomaticResult<Vec<Option<Rank>>> {
        let mut out: Vec<Option<Rank>> = vec![None; self.node_labels.len()];
        // Post-order puts every child before its parent.
        for node in self.newick.ids() {
            let mut lowest = self
                .newick
                .children(node)
                .iter()
                .filter_map(|child| out.get(child.index()).copied().flatten())
                .max();
            if let Some(taxonomy) = &self.label(node)?.taxonomy {
                for taxon in taxonomy.iter_ranks().filter(|taxon| taxon.is_valid_name()) {
                    let rank = taxon.rank()?;
                    if rank != Rank::Species {
                        lowest = lowest.max(Some(rank));
                    }
                }
            }
            if let Some(slot) = out.get_mut(node.index()) {
                *slot = lowest;
            }
        }
        Ok(out)
    }

    /// The highest rank of any named taxon (see [`Taxon::is_named`]) on a node strictly below each
    /// node (leaves included, the node itself excluded), indexed by node. `None` if there is none.
    fn highest_named_rank_below(&self) -> CuratomaticResult<Vec<Option<Rank>>> {
        // The highest named rank on each node or below it.
        let mut at_or_below: Vec<Option<Rank>> = vec![None; self.node_labels.len()];
        let mut out: Vec<Option<Rank>> = vec![None; self.node_labels.len()];
        // Post-order puts every child before its parent.
        for node in self.newick.ids() {
            let below = self
                .newick
                .children(node)
                .iter()
                .filter_map(|child| at_or_below.get(child.index()).copied().flatten())
                .min();
            let mut highest = below;
            if let Some(taxonomy) = &self.label(node)?.taxonomy {
                for taxon in taxonomy.iter_ranks().filter(|taxon| taxon.is_named()) {
                    let rank = taxon.rank()?;
                    highest = Some(highest.map_or(rank, |h| h.min(rank)));
                }
            }
            if let (Some(slot), Some(out_slot)) =
                (at_or_below.get_mut(node.index()), out.get_mut(node.index()))
            {
                *slot = highest;
                *out_slot = below;
            }
        }
        Ok(out)
    }

    /// Give each leaf the taxa that the tree can't hold, from the taxonomy file.
    ///
    /// The file restores the labels of the dummy nodes of an old style tree: each row is the full
    /// taxonomy of a genome, which is the taxa on its path in the tree plus the taxa that were on
    /// its dummy node. This is matched by taxon, never by rank:
    ///
    /// 1. From each leaf, go up the tree to the nearest node with a named taxon (the anchor). Its
    ///    taxa must be in the genome's row, which checks that the tree and the file agree.
    /// 2. Every named taxon of the row that is not on the leaf's path is given to the leaf, e.g. a
    ///    single genome genus (and its species) nested inside a genus on the tree.
    ///
    /// The row is in rank order, not path order, so its position in the row says nothing: a
    /// single genome order nested inside a family on the tree is written before that family.
    /// Removing the dummy nodes and assigning the file gives back exactly the taxa they held.
    ///
    /// The taxa already on the tree are never changed, and placeholders (e.g. `g__`) are never
    /// assigned. A leaf with no named taxon on its path is given its whole row.
    ///
    /// The tree and the file are inputs and are never second guessed: a leaf whose genome is not
    /// in the file, an anchor taxon that is not in the row, and a taxon to assign that is already
    /// on the tree or that more than one leaf would get, are all errors. Reloads the tree.
    pub(crate) fn assign_leaf_taxonomy(
        &mut self,
        file: &TaxonomyFile,
    ) -> CuratomaticResult<LeafAssignment> {
        let conflict = |gid: &GenomeId, message: String| CuratomaticError::TaxonomyFileConflict {
            gid: gid.clone(),
            message,
        };
        let mut out = LeafAssignment::default();
        // The taxa that each leaf gets, and the genome that gets each taxon.
        let mut wanted: Vec<(NodeId, Vec<Taxon>)> = Vec::new();
        let mut wanted_by: HashMap<&Taxon, &GenomeId> = HashMap::new();
        for leaf in self.newick.leaves() {
            let Some(gid) = self.leaf_gid(leaf) else {
                continue;
            };
            let row = file
                .get(gid)
                .ok_or_else(|| CuratomaticError::MissingFromTaxonomyFile {
                    n: 1,
                    first: gid.clone(),
                })?;
            let in_row: BTreeSet<&Taxon> = row.iter_ranks().collect();

            // The named taxa on the path (leaf excluded), checking those of the anchor.
            let mut on_path: BTreeSet<&Taxon> = BTreeSet::new();
            let mut current = self.newick.parent(leaf);
            while let Some(node) = current {
                if let Some(taxonomy) = &self.label(node)?.taxonomy {
                    let is_anchor = on_path.is_empty();
                    for taxon in taxonomy.iter_ranks().filter(|t| t.is_named()) {
                        if is_anchor && !in_row.contains(taxon) {
                            return Err(conflict(
                                gid,
                                format!(
                                    "the tree places it in {taxon}, which is not in its row of \
                                     the taxonomy file"
                                ),
                            ));
                        }
                        on_path.insert(taxon);
                    }
                }
                current = self.newick.parent(node);
            }

            let mut taxa = Vec::new();
            for taxon in row
                .iter_ranks()
                .filter(|t| t.is_named() && !on_path.contains(*t))
            {
                if self.taxon_labels.contains(taxon) {
                    return Err(conflict(
                        gid,
                        format!("{taxon} is for this genome only, but it is already on the tree"),
                    ));
                }
                if let Some(other) = wanted_by.insert(taxon, gid) {
                    return Err(conflict(
                        gid,
                        format!("{taxon} would be given to both {other} and {gid}"),
                    ));
                }
                taxa.push(taxon.clone());
            }
            if !taxa.is_empty() {
                wanted.push((leaf, taxa));
            }
        }

        for (leaf, taxa) in wanted {
            let leaf_taxonomy = &mut self
                .leaves
                .get_mut(&leaf)
                .ok_or_else(|| CuratomaticError::Invariant(format!("node {leaf} is not a leaf")))?
                .taxonomy;
            for taxon in taxa {
                let n = out.n_taxa.entry(taxon.rank()?).or_default();
                *n = n.saturating_add(1);
                leaf_taxonomy.add_taxon(taxon)?;
            }
            out.n_leaves = out.n_leaves.saturating_add(1);
        }

        self.reload()?;
        Ok(out)
    }

    /// The full taxonomy of each leaf, in tree order: every named taxon on the leaf and its
    /// ancestors, in rank order, and in path order (nearest the root first) within a rank, with
    /// the placeholder of each rank that has none (e.g. `g__`). Unlike
    /// [`Tree::leaf_taxonomies`], nested taxa are all kept, so this is the format of the taxonomy
    /// file (`--taxonomy`). Call [`Tree::reload`] first if the tree was changed.
    pub(crate) fn full_taxonomies(&self) -> CuratomaticResult<Vec<(GenomeId, Vec<Taxon>)>> {
        let mut out = Vec::new();
        for leaf in self.newick.leaves() {
            let Some(gid) = self.leaf_gid(leaf) else {
                continue;
            };
            // Read from the leaf up, so each node's taxa go in front of those below it.
            let mut path: Vec<&Taxon> = Vec::new();
            let mut current = Some(leaf);
            while let Some(node) = current {
                if let Some(taxonomy) = &self.label(node)?.taxonomy {
                    let mut taxa: Vec<&Taxon> =
                        taxonomy.iter_ranks().filter(|t| t.is_named()).collect();
                    taxa.append(&mut path);
                    path = taxa;
                }
                current = self.newick.parent(node);
            }
            let mut taxa = Vec::new();
            for rank in Rank::iter() {
                let before = taxa.len();
                for &taxon in &path {
                    if taxon.rank()? == rank {
                        taxa.push(taxon.clone());
                    }
                }
                if taxa.len() == before {
                    taxa.push(Taxon::new(format!("{}__", rank.prefix())));
                }
            }
            out.push((gid.clone(), taxa));
        }
        Ok(out)
    }

    /// The tree in the old style, with a dummy node above each leaf: leaf `G` is written as
    /// `('G|RED=1.0':0.0,'D-G|RED=1.0':0.0)'<the leaf's taxa>|RED=1.0':<its branch length>`, with
    /// the taxa joined by `; `. Every other node is written as it is now. Leaf labels are not
    /// changed: this only formats the tree. Call [`Tree::reload`] first if the tree was changed.
    pub(crate) fn to_zombie_newick(&self) -> String {
        self.newick.to_newick_with(|out, leaf| {
            let Some(info) = self.leaves.get(&leaf) else {
                return;
            };
            let gid = info.gid.as_str();
            let taxa: Vec<&str> = info.taxonomy.iter_ranks().map(Taxon::as_str).collect();
            out.push('(');
            out.push_str(&escape_label(&format!("{gid}|RED=1.0")));
            out.push_str(":0.0,");
            out.push_str(&escape_label(&format!("D-{gid}|RED=1.0")));
            out.push_str(":0.0)");
            out.push_str(&escape_label(&format!("{}|RED=1.0", taxa.join("; "))));
            if let Some(length) = self.newick.length(leaf) {
                out.push(':');
                out.push_str(length);
            }
        })
    }

    /// Check that every leaf holds a taxon, and that one of them is its named species (see the
    /// rules in `AGENTS.md`).
    ///
    /// The species is always the leaf's own: it is never on an internal node, and never a
    /// placeholder.
    pub(crate) fn check_leaves_have_species(&self) -> CuratomaticResult<()> {
        let mut missing = self
            .leaves
            .values()
            .filter(|leaf| {
                !leaf
                    .taxonomy
                    .get_rank(Rank::Species)
                    .iter()
                    .any(Taxon::is_named)
            })
            .map(|leaf| &leaf.gid);
        let Some(first) = missing.next() else {
            return Ok(());
        };
        Err(CuratomaticError::LeafWithoutSpecies {
            n: missing.count().saturating_add(1),
            first: first.clone(),
        })
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

    /// Remove the RED value from the label of every internal node. Leaf labels are never changed.
    pub(crate) fn remove_red_labels(&mut self) -> CuratomaticResult<()> {
        let internal: Vec<NodeId> = self
            .newick
            .ids()
            .filter(|&node| !self.newick.is_leaf(node))
            .collect();
        for node in internal {
            let label = self
                .newick
                .label(node)
                .ok_or_else(|| CuratomaticError::InvalidNodeLabel(format!("node {node}")))?;
            let label = strip_red(label)?.to_owned();
            self.newick.set_internal_label(node, label)?;
        }
        Ok(())
    }

    /// Every taxon that is on the tree, read from the labels of the internal nodes and the
    /// leaves as they are now.
    pub(crate) fn extract_all_taxa_from_tree(&self) -> CuratomaticResult<BTreeSet<Taxon>> {
        let mut out = BTreeSet::new();
        for node in self.newick.ids() {
            out.extend(self.node_taxonomy(node)?.iter_ranks().cloned());
        }
        Ok(out)
    }

    /// The number of nodes that hold each named taxon, read from the labels of the internal
    /// nodes and the leaves as they are now.
    pub(crate) fn taxon_counts(&self) -> CuratomaticResult<TaxonCounts> {
        let mut taxa = Vec::new();
        for node in self.newick.ids() {
            taxa.extend(self.node_taxonomy(node)?.iter_ranks().cloned());
        }
        Ok(TaxonCounts::from_taxa(&taxa))
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
        let species_genera = self.species_genera_of_new_genomes(rank, &gid_to_taxonomy);

        let child_rank = rank.child().ok_or_else(|| {
            CuratomaticError::Invariant("the species cannot be refined, it has no child".to_owned())
        })?;
        let median_red = red_dict.get(rank)?;
        let parent_red_thresh = red_dict.get(parent_rank(rank)?)?.midpoint(median_red);
        let child_red_thresh = red_dict.get(child_rank)?.midpoint(median_red);
        let highest_below = self.highest_named_rank_below()?;

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
                            // Stop if a taxon of a higher rank is below this node (the original
                            // node included): the taxon would be placed above it. Every higher
                            // ancestor has it below too, so the walk ends here.
                            if highest_below
                                .get(cur_node.index())
                                .copied()
                                .flatten()
                                .is_some_and(|higher| higher < rank)
                            {
                                break;
                            }

                            if median_red.distance(cur_red) < original_red_delta {
                                // This brings the RED closer to the median.
                                out.push(taxon, cur_node, node);
                            } else if original_bs < min_bs {
                                // This moves it further from the median, only do so to
                                // improve on a node that is not supported: support is traded
                                // for RED, while staying within the RED range of the rank.
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
                let child_taxa =
                    self.child_rank_taxa(node, rank, &gid_to_taxonomy, &species_genera)?;

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

                    let child_taxa =
                        self.child_rank_taxa(cur_node, rank, &gid_to_taxonomy, &species_genera)?;
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

    /// For the genus only: the genus that the species name of each new genome gives it, keyed by
    /// genome ID. A new genome is one without a prior taxonomy (it was not in the previous
    /// release), and its species is held by its leaf (from `--taxonomy`), e.g. `g__Noviluteimonas`
    /// for `s__Noviluteimonas sp057139745` (see [`Taxon::genus_of_species`]).
    ///
    /// Moving a genus down follows the genomes that the prior taxonomy puts in it. A new genome is
    /// not in the prior, but a binomial species name says which genus it belongs to, so the genus
    /// is not moved down away from it: a genus is never split from the genomes its species name.
    fn species_genera_of_new_genomes(
        &self,
        rank: Rank,
        gid_to_taxonomy: &BTreeMap<&GenomeId, &Taxonomy>,
    ) -> BTreeMap<GenomeId, Taxon> {
        if rank != Rank::Genus {
            return BTreeMap::new();
        }
        self.leaves
            .values()
            .filter(|leaf| !gid_to_taxonomy.contains_key(&leaf.gid))
            .filter_map(|leaf| {
                let genus = leaf
                    .taxonomy
                    .get_rank(Rank::Species)
                    .iter()
                    .find_map(Taxon::genus_of_species)?;
                Some((leaf.gid.clone(), genus))
            })
            .collect()
    }

    /// The taxa of `rank` that the metadata assigns to the genomes at or below each child of
    /// `node` (a child may be a leaf itself).
    ///
    /// A new genome (one without a prior taxonomy) is put in the genus that its species name
    /// gives, see [`Tree::species_genera_of_new_genomes`].
    fn child_rank_taxa<'a>(
        &self,
        node: NodeId,
        rank: Rank,
        gid_to_taxonomy: &BTreeMap<&GenomeId, &'a Taxonomy>,
        species_genera: &'a BTreeMap<GenomeId, Taxon>,
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
                } else if let Some(genus) = gid.and_then(|gid| species_genera.get(gid)) {
                    out.entry(child).or_default().insert(genus);
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
            transfer_taxon(self, original_node, best_node, taxon)?;
            moves.push(Move {
                taxon: taxon.clone(),
                from: original_node,
                to: best_node,
            });
        }
        Ok(moves)
    }
}

/// Move a taxon from one node to another. Either node may be a leaf, whose label is never changed:
/// its taxa are on its [`Leaf`].
fn transfer_taxon(
    tree: &mut Tree,
    old_node: NodeId,
    new_node: NodeId,
    taxon: &Taxon,
) -> CuratomaticResult<()> {
    let mut old_taxonomy = tree.node_taxonomy(old_node)?;
    old_taxonomy.remove_taxon(taxon)?;

    // An internal node without a bootstrap is labelled by its taxonomy alone, so removing its last
    // taxon would leave nothing to label it with (a leaf keeps its genome ID). Neither node has
    // been changed at this point.
    if !tree.leaves.contains_key(&old_node)
        && NodeLabel::from_node(&tree.newick, old_node)?
            .bootstrap
            .is_none()
        && old_taxonomy.iter_ranks().next().is_none()
    {
        return Err(CuratomaticError::Invariant(format!(
            "moving {taxon} would leave node {old_node} with neither a bootstrap nor a taxonomy"
        )));
    }

    tree.set_node_taxonomy(old_node, old_taxonomy)?;
    tree.add_node_taxon(new_node, taxon.clone())
}

/// `None` for a taxonomy without any taxa.
fn non_empty(taxonomy: Taxonomy) -> Option<Taxonomy> {
    let has_taxa = taxonomy.iter_ranks().next().is_some();
    has_taxa.then_some(taxonomy)
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
/// A novel taxon is only created when no valid name (see [`Taxon::is_valid_name`]) of a lower
/// rank, other than a species, is on the highest node or below it. Otherwise its leaves keep the
/// placeholder for that rank only: each lower missing rank is decided on its own candidate node.
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

    // Read before any novel taxon is created, so only the names already on the tree count.
    let lowest_valid_below = tree.lowest_valid_rank_below()?;

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
        let mut n_skipped: usize = 0;
        for (candidate, leaves) in &candidate_to_leaves {
            // A curator never places a novel taxon above a valid name of a lower rank. Only this
            // rank is skipped: each lower rank is decided on its own candidate node.
            let lowest = lowest_valid_below.get(candidate.index()).copied().flatten();
            if lowest.is_some_and(|lowest| lowest > *rank) {
                n_skipped = n_skipped.saturating_add(1);
                continue;
            }

            let gids = leaves
                .iter()
                .filter_map(|&leaf| tree.leaf_gid(leaf).cloned())
                .collect::<BTreeSet<GenomeId>>();
            let new_taxon = generate_unique_rank(*rank, &gids, metadata, &known_taxa)?;
            known_taxa.insert(new_taxon.clone());
            new_taxa_created.insert(new_taxon.clone());

            // Add the taxon to the node (the leaf itself may be the candidate).
            tree.add_node_taxon(*candidate, new_taxon)?;
        }
        info!(
            "Skipped {n_skipped} novel {} taxa that have a valid name below them.",
            rank.prefix()
        );
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
