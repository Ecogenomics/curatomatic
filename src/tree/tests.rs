use std::path::Path;

use super::*;

fn tree(newick: &str) -> CuratomaticResult<Tree> {
    Tree::from_newick(NewickTree::parse(newick)?)
}

/// The node at a position in post-order.
fn node(tree: &Tree, index: usize) -> NodeId {
    tree.newick.ids().nth(index).unwrap()
}

fn label(tree: &Tree, index: usize) -> &str {
    tree.newick.label(node(tree, index)).unwrap()
}

/// The taxa held by the leaf at a position in post-order, whose label must be unchanged.
fn leaf_taxa(tree: &Tree, index: usize) -> String {
    let leaf = &tree.leaves[&node(tree, index)];
    assert_eq!(
        label(tree, index).split('|').next(),
        Some(leaf.gid.as_str())
    );
    leaf.taxonomy.to_string()
}

/// A tree whose leaves hold the given taxa, which a leaf label never can.
fn tree_with_leaf_taxa(newick: &str, leaf_taxa: &[(usize, &str)]) -> CuratomaticResult<Tree> {
    let mut tree = tree(newick)?;
    for &(index, taxa) in leaf_taxa {
        let leaf = tree.leaves.get_mut(&node(&tree, index)).unwrap();
        leaf.taxonomy = Taxonomy::from_string(taxa)?;
    }
    tree.reload()?;
    Ok(tree)
}

fn red_dict() -> Red {
    Red::from_tsv("rank\tmedian\nPhylum\t0.3\nClass\t0.4\nOrder\t0.5\nFamily\t0.7\nGenus\t0.9\n")
        .unwrap()
}

/// Post-order: G1 (0), G3 (1), their species node (2), G2 (3), its species node (4), root (5).
const TWO_GENOMES: &str = "((G1|RED=1.000:0.1,G3|RED=1.000:0.1)'g__A; s__A one|RED=1.000':0.1,\
     (G2|RED=1.000:0.1)'s__B b|RED=1.000':0.1)'90.0:f__F|RED=0.700';";

#[test]
fn reading_a_tree_does_not_rewrite_its_labels() {
    let tree = tree(TWO_GENOMES).unwrap();
    assert_eq!(label(&tree, 2), "g__A; s__A one|RED=1.000");
    assert_eq!(label(&tree, 4), "s__B b|RED=1.000");
    assert_eq!(label(&tree, 5), "90.0:f__F|RED=0.700");
}

#[test]
fn transfer_taxon_moves_the_taxon_and_rewrites_both_labels() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, species) = (node(&tree, 5), node(&tree, 4));
    transfer_taxon(&mut tree, root, species, &Taxon::new("f__F")).unwrap();
    assert_eq!(label(&tree, 5), "90.0|RED=0.700");
    assert_eq!(label(&tree, 4), "f__F; s__B b|RED=1.000");
}

#[test]
fn transfer_taxon_does_not_strip_a_node_that_has_no_bootstrap() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, species) = (node(&tree, 5), node(&tree, 4));
    assert!(matches!(
        transfer_taxon(&mut tree, species, root, &Taxon::new("s__B b")),
        Err(CuratomaticError::Invariant(_))
    ));
    // Nothing was changed.
    assert_eq!(label(&tree, 4), "s__B b|RED=1.000");
    assert_eq!(label(&tree, 5), "90.0:f__F|RED=0.700");
}

#[test]
fn transfer_taxon_keeps_a_node_without_a_bootstrap_that_has_another_taxon() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, species) = (node(&tree, 5), node(&tree, 2));
    transfer_taxon(&mut tree, species, root, &Taxon::new("g__A")).unwrap();
    assert_eq!(label(&tree, 2), "s__A one|RED=1.000");
    assert_eq!(label(&tree, 5), "90.0:f__F; g__A|RED=0.700");
}

#[test]
fn transfer_taxon_moves_a_taxon_to_and_from_a_leaf() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, leaf) = (node(&tree, 5), node(&tree, 3));
    transfer_taxon(&mut tree, root, leaf, &Taxon::new("f__F")).unwrap();
    assert_eq!(label(&tree, 5), "90.0|RED=0.700");
    assert_eq!(leaf_taxa(&tree, 3), "f__F");

    // A leaf keeps its genome ID, so its last taxon can be moved away.
    transfer_taxon(&mut tree, leaf, root, &Taxon::new("f__F")).unwrap();
    assert_eq!(label(&tree, 3), "G2|RED=1.000");
    assert_eq!(label(&tree, 5), "90.0:f__F|RED=0.700");
}

#[test]
fn a_taxon_on_a_leaf_is_held_by_the_leaf_and_never_by_its_label() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, leaf) = (node(&tree, 5), node(&tree, 3));
    transfer_taxon(&mut tree, root, leaf, &Taxon::new("f__F")).unwrap();
    tree.reload().unwrap();
    assert_eq!(label(&tree, 3), "G2|RED=1.000");
    assert_eq!(
        tree.label(leaf).unwrap().taxonomy,
        Some(Taxonomy::from_string("f__F").unwrap())
    );
    assert!(tree.taxon_nodes().contains(&(Taxon::new("f__F"), leaf)));
}

#[test]
fn removing_the_red_labels_never_changes_a_leaf_label() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    tree.remove_red_labels().unwrap();
    let labels: Vec<_> = (0..6).map(|i| label(&tree, i)).collect();
    assert_eq!(
        labels,
        [
            "G1|RED=1.000",
            "G3|RED=1.000",
            "g__A; s__A one",
            "G2|RED=1.000",
            "s__B b",
            "90.0:f__F"
        ]
    );
}

#[test]
fn a_novel_domain_is_an_error_as_the_domain_has_no_parent_rank() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let missing = tree.get_leaf_nodes_missing_rank().unwrap();
    assert!(missing.values().all(|ranks| ranks.contains(&Rank::Domain)));

    let known_taxa = tree.taxon_labels.clone();
    assert_eq!(
        fill_in_missing_ranks(
            &mut tree,
            &missing,
            BootstrapValue::new(95.0),
            &red_dict(),
            &Metadata::default(),
            known_taxa,
        ),
        Err(CuratomaticError::NoParentRank('d'))
    );
}

#[test]
fn refining_the_domain_is_an_error_as_it_has_no_parent_rank() {
    let tree = tree(TWO_GENOMES).unwrap();
    assert!(matches!(
        tree.refine_existing_labels(
            Rank::Domain,
            BootstrapValue::new(95.0),
            &red_dict(),
            &Metadata::default()
        ),
        Err(CuratomaticError::NoParentRank('d'))
    ));
}

/// `g__B` (bootstrap 50) sits below an unnamed node Y (RED 0.88, bootstrap 100), which is closer to
/// the genus median and would otherwise qualify. Above Y is X, labelled `{x}`, with G4 below it.
fn nested_tree(x: &str) -> String {
    format!(
        "(((((G1|RED=1.000:0.1)'s__B one|RED=1.000':0.1,(G2|RED=1.000:0.1)'s__B two|RED=1.000':0.1)\
         '50.0:g__B|RED=0.930':0.1,(G3|RED=1.000:0.1)'s__A three|RED=1.000':0.1)'100.0|RED=0.880':0.1,\
         (G4|RED=1.000:0.1)'s__A four|RED=1.000':0.1)'{x}|RED=0.820':0.1,\
         (G5|RED=1.000:0.1)'s__C five|RED=1.000':0.1)'100.0:f__F|RED=0.700';"
    )
}

/// `g__B` (bootstrap 50) sits directly below the node of `g__A` (RED 0.88, bootstrap 100), which is
/// closer to the genus median and would otherwise qualify.
const NESTED_ON_ANCESTOR_NODE: &str = "((((G1|RED=1.000:0.1)'s__B one|RED=1.000':0.1,\
     (G2|RED=1.000:0.1)'s__B two|RED=1.000':0.1)'50.0:g__B|RED=0.930':0.1,\
     (G4|RED=1.000:0.1)'s__A four|RED=1.000':0.1)'100.0:g__A|RED=0.880':0.1,\
     (G5|RED=1.000:0.1)'s__C five|RED=1.000':0.1)'100.0:f__F|RED=0.700';";

/// Refine the genera of a tree, as is done for each rank when curatomatic runs.
fn refine_genus(newick: &str) -> Tree {
    let mut tree = tree(newick).unwrap();
    let min_bs = BootstrapValue::new(95.0);
    let refinement = tree
        .refine_existing_labels(Rank::Genus, min_bs, &red_dict(), &Metadata::default())
        .unwrap();
    tree.process_refinement(&refinement, Rank::Genus, min_bs, &red_dict())
        .unwrap();
    tree.reload().unwrap();
    tree
}

/// The label of the only node whose label contains `needle`.
fn label_with<'a>(tree: &'a Tree, needle: &str) -> &'a str {
    let mut found = tree
        .newick
        .ids()
        .filter_map(|node| tree.newick.label(node))
        .filter(|l| l.contains(needle));
    let label = found.next().unwrap();
    assert!(
        found.next().is_none(),
        "more than one label contains {needle}"
    );
    label
}

#[test]
fn a_taxon_that_is_not_nested_is_moved_up_to_the_qualifying_node() {
    // Without g__A the fixture qualifies, so the tests below show that nesting is what stops it.
    let tree = refine_genus(&nested_tree("100.0"));
    assert_eq!(label_with(&tree, "g__B"), "100.0:g__B|RED=0.880");
    assert_eq!(label_with(&tree, "RED=0.930"), "50.0|RED=0.930");
}

#[test]
fn a_nested_taxon_is_not_moved_up_part_of_the_way_to_its_ancestor() {
    let tree = refine_genus(&nested_tree("100.0:g__A"));
    assert_eq!(label_with(&tree, "g__B"), "50.0:g__B|RED=0.930");
    assert_eq!(label_with(&tree, "g__A"), "100.0:g__A|RED=0.820");
    assert_eq!(label_with(&tree, "RED=0.880"), "100.0|RED=0.880");
}

#[test]
fn a_nested_taxon_is_not_moved_onto_the_node_of_its_ancestor() {
    let tree = refine_genus(NESTED_ON_ANCESTOR_NODE);
    assert_eq!(label_with(&tree, "g__B"), "50.0:g__B|RED=0.930");
    assert_eq!(label_with(&tree, "g__A"), "100.0:g__A|RED=0.880");
}

#[test]
fn nested_inside_finds_the_closest_named_ancestor_of_the_same_rank() {
    let tree = tree(&nested_tree("100.0:g__A; g__")).unwrap();
    let b = tree
        .newick
        .ids()
        .find(|&node| tree.newick.label(node).is_some_and(|l| l.contains("g__B")))
        .unwrap();
    let (genus, family) = (Taxon::new("g__B"), Taxon::new("f__F"));
    assert_eq!(
        tree.nested_inside(b, Rank::Genus, &genus),
        Some(&Taxon::new("g__A"))
    );
    assert_eq!(tree.nested_inside(b, Rank::Family, &family), None);
}

// ---- The taxonomy file (`--taxonomy`) ----

fn taxonomy_file(contents: &str) -> TaxonomyFile {
    TaxonomyFile::from_reader(contents.as_bytes(), Path::new("taxonomy.tsv")).unwrap()
}

/// A tree in the r237 format: leaves without a RED, and `|red=`. Post-order: G1 (0), G2 (1),
/// `g__A` (2), G3 (3), `f__F` (4), G4 (5), `c__C` (6), root (7).
const NEW_STYLE: &str = "(((G1:0.1,G2:0.1)'100:g__A|red=0.95':0.1,G3:0.1)'100:f__F|red=0.70':0.1,\
     (G4:0.1)'100:c__C|red=0.40':0.1)'d__D|red=0.0000';";

/// The taxonomy file for [`NEW_STYLE`], where the file disagrees with the tree for G1's genus,
/// G3's genus is a placeholder, and G4 has a phylum that the tree doesn't.
const NEW_STYLE_TAXONOMY: &str = "\
G1\td__D;p__P;c__C;o__O;f__F;g__Other;s__A one
G2\td__D;p__P;c__C;o__O;f__F;g__A;s__{unresolved} sp2
G3\td__D;p__P;c__C;o__O;f__F;g__;s__B three
G4\td__D;p__Q;c__C;o__O2;f__F2;g__G;s__G four
G9\td__D;p__P;c__C;o__O;f__F;g__A;s__A nine
";

fn assigned(newick: &str, taxonomy: &str) -> (Tree, LeafAssignment) {
    let mut tree = tree(newick).unwrap();
    let assignment = tree.assign_leaf_taxonomy(&taxonomy_file(taxonomy)).unwrap();
    (tree, assignment)
}

#[test]
fn a_leaf_gets_the_ranks_below_the_lowest_rank_on_its_path() {
    let (tree, assignment) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    // The lowest rank above G1 and G2 is the genus: only the species is assigned.
    assert_eq!(leaf_taxa(&tree, 0), "s__A one");
    // The lowest rank above G4 is the class: the order down to the species.
    assert_eq!(leaf_taxa(&tree, 5), "o__O2; f__F2; g__G; s__G four");
    assert_eq!(assignment.n_leaves, 4);
    assert_eq!(
        assignment.n_taxa,
        BTreeMap::from([
            (Rank::Order, 1),
            (Rank::Family, 1),
            (Rank::Genus, 1),
            (Rank::Species, 4)
        ])
    );
    assert!(assignment.skipped.is_empty());
}

#[test]
fn the_taxa_on_the_tree_are_never_changed_by_the_taxonomy_file() {
    let before = tree(NEW_STYLE).unwrap();
    let (after, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    // G1 is in `g__Other` in the file, the tree's `g__A` wins.
    assert!(!after.taxon_labels.contains(&Taxon::new("g__Other")));
    for index in [2, 4, 6, 7] {
        assert_eq!(label(&after, index), label(&before, index));
    }
    assert!(after.taxon_labels.is_superset(&before.taxon_labels));
}

#[test]
fn a_placeholder_in_the_taxonomy_file_is_never_assigned() {
    let (tree, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    // G3's genus is `g__` in the file.
    assert_eq!(leaf_taxa(&tree, 3), "s__B three");
}

#[test]
fn a_gap_above_the_lowest_rank_is_not_filled() {
    let (tree, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    // G4 has no phylum on its path, but its class is the lowest rank, so `p__Q` is not assigned.
    assert!(!tree.taxon_labels.contains(&Taxon::new("p__Q")));
    let missing = tree.get_leaf_nodes_missing_rank().unwrap();
    assert_eq!(
        missing.get(&node(&tree, 5)),
        Some(&BTreeSet::from([Rank::Phylum]))
    );
}

#[test]
fn an_unresolved_species_reaches_the_leaf_unchanged() {
    let (tree, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    assert_eq!(leaf_taxa(&tree, 1), "s__{unresolved} sp2");
    let taxonomies = tree.leaf_taxonomies().unwrap();
    assert_eq!(
        taxonomies[1].1.last(),
        Some(&Taxon::new("s__{unresolved} sp2"))
    );
}

#[test]
fn an_old_style_tree_gets_nothing_from_the_taxonomy_file() {
    let taxonomy = "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A other\n\
                    G2\td__D;p__P;c__C;o__O;f__F;g__B;s__B other\n\
                    G3\td__D;p__P;c__C;o__O;f__F;g__A;s__A third\n";
    let (tree, assignment) = assigned(TWO_GENOMES, taxonomy);
    assert_eq!(assignment, LeafAssignment::default());
    assert_eq!(leaf_taxa(&tree, 0), "");
    assert_eq!(leaf_taxa(&tree, 3), "");
}

#[test]
fn a_leaf_missing_from_the_taxonomy_file_is_an_error() {
    let mut tree = tree(NEW_STYLE).unwrap();
    let file = taxonomy_file("G1\td__D;s__A one\n");
    assert_eq!(
        tree.assign_leaf_taxonomy(&file),
        Err(CuratomaticError::MissingFromTaxonomyFile {
            n: 1,
            first: GenomeId::new("G2")
        })
    );
}

#[test]
fn a_taxon_already_on_the_tree_or_wanted_by_two_leaves_is_skipped() {
    // G3 and G4 both want `g__G`, and G4 wants `f__F`, which is already on the tree.
    let taxonomy = "\
G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one
G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two
G3\td__D;p__P;c__C;o__O;f__F;g__G;s__G three
G4\td__D;p__P;c__C;o__O2;f__F;g__G;s__G four
";
    let (tree, assignment) = assigned(NEW_STYLE, taxonomy);
    assert_eq!(
        assignment.skipped,
        BTreeSet::from([Taxon::new("f__F"), Taxon::new("g__G")])
    );
    assert_eq!(leaf_taxa(&tree, 3), "s__G three");
    assert_eq!(leaf_taxa(&tree, 5), "o__O2; s__G four");
}

#[test]
fn leaf_taxonomies_read_the_ranks_in_order_from_the_leaf_and_its_ancestors() {
    let (tree, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    let taxonomies = tree.leaf_taxonomies().unwrap();
    let rows: Vec<(&str, String)> = taxonomies
        .iter()
        .map(|(gid, taxa)| {
            let taxa: Vec<&str> = taxa.iter().map(Taxon::as_str).collect();
            (gid.as_str(), taxa.join(";"))
        })
        .collect();
    // A rank that no node holds is a placeholder.
    assert_eq!(
        rows,
        vec![
            ("G1", "d__D;p__;c__;o__;f__F;g__A;s__A one".to_owned()),
            (
                "G2",
                "d__D;p__;c__;o__;f__F;g__A;s__{unresolved} sp2".to_owned()
            ),
            ("G3", "d__D;p__;c__;o__;f__F;g__;s__B three".to_owned()),
            ("G4", "d__D;p__;c__C;o__O2;f__F2;g__G;s__G four".to_owned()),
        ]
    );
}

#[test]
fn leaf_taxonomies_use_the_nearest_of_two_nested_taxa() {
    let tree = tree_with_leaf_taxa(
        "(((G1:0.1,G2:0.1)'100:g__Inner|red=0.95':0.1,G3:0.1)'100:g__Outer|red=0.9':0.1)\
         'd__D; g__|red=0.0';",
        &[(0, "s__A one"), (1, "s__A two"), (3, "s__B three")],
    )
    .unwrap();
    let taxonomies = tree.leaf_taxonomies().unwrap();
    let genus = |idx: usize| taxonomies[idx].1[5].as_str();
    assert_eq!(genus(0), "g__Inner");
    assert_eq!(genus(2), "g__Outer");
}

#[test]
fn a_genus_is_moved_down_onto_a_leaf() {
    // `g__A` (bootstrap 50) is on a node above G1 and G2, but only G1 is in it.
    let mut tree =
        tree("((G1:0.1,G2:0.1)'50:g__A|red=0.95':0.1,G3:0.1)'100:f__F|red=0.70';").unwrap();
    let mut metadata = Metadata::default();
    for (gid, taxonomy) in [
        ("G1", "f__F;g__A"),
        ("G2", "f__F;g__B"),
        ("G3", "f__F;g__C"),
    ] {
        metadata.set_taxonomy(GenomeId::new(gid), Taxonomy::from_string(taxonomy).unwrap());
    }
    let min_bs = BootstrapValue::new(95.0);
    let refinement = tree
        .refine_existing_labels(Rank::Genus, min_bs, &red_dict(), &metadata)
        .unwrap();
    let moves = tree
        .process_refinement(&refinement, Rank::Genus, min_bs, &red_dict())
        .unwrap();
    tree.reload().unwrap();
    assert_eq!(moves.len(), 1);
    assert_eq!(leaf_taxa(&tree, 0), "g__A");
    assert_eq!(label(&tree, 2), "50.0|RED=0.950");
}

#[test]
fn a_novel_taxon_can_be_placed_on_a_leaf() {
    // Neither node above G1 is supported, so the novel genus can only go on the leaf.
    let mut tree = tree(
        "(((G1:0.1,G2:0.1)'50:red=0.97':0.1,G3:0.1)'50:f__F|red=0.70':0.1,\
         (G4:0.1)'100:f__X;g__X|red=0.5')'d__D;p__P;c__C;o__O|red=0.0';",
    )
    .unwrap();
    let file = taxonomy_file(
        "G1\td__D;s__{unresolved} G1\nG2\td__D;s__{unresolved} G2\n\
         G3\td__D;s__{unresolved} G3\nG4\td__D;s__X four\n",
    );
    tree.assign_leaf_taxonomy(&file).unwrap();
    let missing = tree.get_leaf_nodes_missing_rank().unwrap();
    let known = tree.taxon_labels.clone();
    let novel = fill_in_missing_ranks(
        &mut tree,
        &missing,
        BootstrapValue::new(95.0),
        &red_dict(),
        &Metadata::default(),
        known,
    )
    .unwrap();
    assert_eq!(
        novel,
        BTreeSet::from([
            Taxon::new("g__G1"),
            Taxon::new("g__G2"),
            Taxon::new("g__G3")
        ])
    );
    assert_eq!(leaf_taxa(&tree, 0), "g__G1; s__{unresolved} G1");
}

#[test]
fn an_unresolvable_taxon_may_be_on_more_than_one_node() {
    let tree = tree_with_leaf_taxa(
        "((G1:0.1,G2:0.1)'100:f__F|red=0.7')'d__D|red=0.0';",
        &[
            (0, "g__UNRESOLVABLE; s__A one"),
            (1, "g__UNRESOLVABLE; s__A two"),
        ],
    )
    .unwrap();
    assert!(
        tree.taxon_labels
            .contains(&Taxon::unresolvable(Rank::Genus))
    );
}

#[test]
fn any_other_taxon_on_more_than_one_node_is_an_error() {
    assert!(matches!(
        tree_with_leaf_taxa(
            "((G1:0.1,G2:0.1)'100:f__F|red=0.7')'d__D|red=0.0';",
            &[(0, "g__X"), (1, "g__X")]
        ),
        Err(CuratomaticError::NonUniqueTaxon(_))
    ));
}

// ---- Novel taxa are only created where nothing below is named ----

/// Fill in the missing ranks of a tree whose leaves hold the given taxa, as `main` does after
/// refinement. Returns the tree and the novel taxa.
fn fill_novel(newick: &str, leaf_taxa: &[(usize, &str)]) -> (Tree, BTreeSet<Taxon>) {
    let mut tree = tree_with_leaf_taxa(newick, leaf_taxa).unwrap();
    let missing = tree.get_leaf_nodes_missing_rank().unwrap();
    let known = tree.taxon_labels.clone();
    let novel = fill_in_missing_ranks(
        &mut tree,
        &missing,
        BootstrapValue::new(95.0),
        &red_dict(),
        &Metadata::default(),
        known,
    )
    .unwrap();
    tree.reload().unwrap();
    (tree, novel)
}

/// The taxonomy written for each genome, joined by `;`.
fn written_taxonomies(tree: &Tree) -> Vec<String> {
    tree.leaf_taxonomies()
        .unwrap()
        .iter()
        .map(|(gid, taxa)| {
            let taxa: Vec<&str> = taxa.iter().map(Taxon::as_str).collect();
            format!("{gid}\t{}", taxa.join(";"))
        })
        .collect()
}

/// G1 (0) and G2 (1) below a supported node X (RED 0.95), below the root, which holds every rank
/// above the genus.
const MISSING_GENUS: &str =
    "((G1:0.1,G2:0.1)'100|red=0.95':0.1)'d__D;p__P;c__C;o__O;f__F|red=0.0';";

#[test]
fn a_missing_phylum_above_a_named_class_is_not_created() {
    // The highest phylum candidate (RED 0.3) is above `c__C`.
    let (tree, novel) = fill_novel(
        "(((G1:0.1,G2:0.1)'100:c__C;o__O;f__F;g__G|red=0.5':0.1)'100|red=0.3':0.1)'d__D|red=0.0';",
        &[(0, "s__G one"), (1, "s__G two")],
    );
    assert!(novel.is_empty());
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__;c__C;o__O;f__F;g__G;s__G one",
            "G2\td__D;p__;c__C;o__O;f__F;g__G;s__G two"
        ]
    );
}

#[test]
fn a_missing_genus_above_only_unresolved_species_is_created() {
    let (tree, novel) = fill_novel(
        MISSING_GENUS,
        &[(0, "s__{unresolved} G1"), (1, "s__{unresolved} G2")],
    );
    assert_eq!(novel, BTreeSet::from([Taxon::new("g__G1")]));
    assert_eq!(label(&tree, 2), "100.0:g__G1|RED=0.950");
}

#[test]
fn a_missing_genus_above_a_named_species_is_not_created() {
    let (tree, novel) = fill_novel(
        MISSING_GENUS,
        &[(0, "s__{unresolved} G1"), (1, "s__Chlorobium limicola")],
    );
    assert!(novel.is_empty());
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__P;c__C;o__O;f__F;g__;s__{unresolved} G1",
            "G2\td__D;p__P;c__C;o__O;f__F;g__;s__Chlorobium limicola"
        ]
    );
}

#[test]
fn a_wgs_style_name_below_blocks_creation() {
    // The highest family candidate (RED 0.75) is above `g__CAYYKD01`, and G3's genus is skipped
    // with its family.
    let (_, novel) = fill_novel(
        "(((G1:0.1,G2:0.1)'100:g__CAYYKD01|red=0.95':0.1,G3:0.1)'100|red=0.75':0.1)\
         'd__D;p__P;c__C;o__O|red=0.0';",
        &[
            (0, "s__{unresolved} G1"),
            (1, "s__{unresolved} G2"),
            (3, "s__{unresolved} G3"),
        ],
    );
    assert!(novel.is_empty());
}

#[test]
fn a_lower_rank_on_the_candidate_node_itself_blocks_creation() {
    // The highest phylum candidate is the node of `c__Chlorobiia`.
    let (tree, novel) = fill_novel(
        "((G1:0.1,G2:0.1)'100:c__Chlorobiia;o__O;f__F;g__G|red=0.4':0.1)'d__D|red=0.0';",
        &[(0, "s__{unresolved} G1"), (1, "s__{unresolved} G2")],
    );
    assert!(novel.is_empty());
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__;c__Chlorobiia;o__O;f__F;g__G;s__{unresolved} G1",
            "G2\td__D;p__;c__Chlorobiia;o__O;f__F;g__G;s__{unresolved} G2"
        ]
    );
}

#[test]
fn the_lower_ranks_of_a_skipped_genome_are_not_filled() {
    // The phylum is skipped at X (RED 0.2), above `c__C`. Nothing named is below Z (RED 0.75), the
    // highest family candidate, but the family and genus are not filled.
    let (tree, novel) = fill_novel(
        "((((G1:0.1,G2:0.1)'100|red=0.75':0.1)'100:c__C;o__O|red=0.45':0.1)'100|red=0.2':0.1)\
         'd__D|red=0.0';",
        &[(0, "s__{unresolved} G1"), (1, "s__{unresolved} G2")],
    );
    assert!(novel.is_empty());
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__;c__C;o__O;f__;g__;s__{unresolved} G1",
            "G2\td__D;p__;c__C;o__O;f__;g__;s__{unresolved} G2"
        ]
    );
}

#[test]
fn a_genome_that_is_not_skipped_still_gets_its_lower_ranks() {
    // Without `c__C; o__O` the same tree gets a novel phylum and family (and more).
    let (_, novel) = fill_novel(
        "((((G1:0.1,G2:0.1)'100|red=0.75':0.1)'100|red=0.45':0.1)'100|red=0.2':0.1)'d__D|red=0.0';",
        &[(0, "s__{unresolved} G1"), (1, "s__{unresolved} G2")],
    );
    assert!(novel.contains(&Taxon::new("p__G1")));
    assert!(novel.contains(&Taxon::new("f__G1")));
}
