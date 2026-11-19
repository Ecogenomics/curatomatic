use super::*;

fn tree(newick: &str) -> CuratomaticResult<Tree> {
    Tree::from_newick(NewickTree::parse(newick)?)
}

/// The node at a position in post-order.
fn node(tree: &Tree, index: usize) -> NodeId {
    tree.tree.ids().nth(index).unwrap()
}

fn label(tree: &Tree, index: usize) -> &str {
    tree.tree.label(node(tree, index)).unwrap()
}

fn red_dict() -> Red {
    Red::from_json(
        r#"{"phylum": 0.3, "class": 0.4, "order": 0.5, "family": 0.7, "genus": 0.9, "species": 1.0}"#,
    )
    .unwrap()
}

/// Post-order: G1 (0), D-G1 (1), their species node (2), G2 (3), its species node (4), root (5).
const TWO_GENOMES: &str = "((G1|RED=1.000:0.1,D-G1|RED=1.000:0.1)'g__A; s__|RED=1.000':0.1,\
     (G2|RED=1.000:0.1)'s__B b|RED=1.000':0.1)'90.0:f__F|RED=0.700';";

#[test]
fn a_species_placeholder_is_replaced_and_the_other_taxa_are_kept() {
    let tree = tree(TWO_GENOMES).unwrap();
    assert_eq!(label(&tree, 2), "g__A; s__{unresolved} G1|RED=1.000");
    assert!(tree.taxon_labels.contains(&Taxon::new("g__A")));
    assert!(tree.taxon_labels.contains(&Taxon::new("s__{unresolved} G1")));
}

#[test]
fn a_species_placeholder_on_its_own_is_replaced() {
    let tree = tree("((G1|RED=1.000:0.1)'s__|RED=1.000':0.1)'90.0:g__A|RED=0.900';").unwrap();
    assert_eq!(label(&tree, 1), "s__{unresolved} G1|RED=1.000");
}

#[test]
fn a_label_without_a_species_placeholder_is_not_rewritten() {
    let tree = tree(TWO_GENOMES).unwrap();
    assert_eq!(label(&tree, 4), "s__B b|RED=1.000");
    assert_eq!(label(&tree, 5), "90.0:f__F|RED=0.700");
}

#[test]
fn transfer_taxon_moves_the_taxon_and_rewrites_both_labels() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, species) = (node(&tree, 5), node(&tree, 4));
    transfer_taxon(&mut tree.tree, root, species, &Taxon::new("f__F")).unwrap();
    assert_eq!(label(&tree, 5), "90.0|RED=0.700");
    assert_eq!(label(&tree, 4), "f__F; s__B b|RED=1.000");
}

#[test]
fn transfer_taxon_does_not_strip_a_node_that_has_no_bootstrap() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, species) = (node(&tree, 5), node(&tree, 4));
    assert!(matches!(
        transfer_taxon(&mut tree.tree, species, root, &Taxon::new("s__B b")),
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
    transfer_taxon(&mut tree.tree, species, root, &Taxon::new("g__A")).unwrap();
    assert_eq!(label(&tree, 2), "s__{unresolved} G1|RED=1.000");
    assert_eq!(label(&tree, 5), "90.0:f__F; g__A|RED=0.700");
}

#[test]
fn transfer_taxon_does_not_move_a_taxon_to_a_leaf() {
    let mut tree = tree(TWO_GENOMES).unwrap();
    let (root, leaf) = (node(&tree, 5), node(&tree, 3));
    assert!(matches!(
        transfer_taxon(&mut tree.tree, root, leaf, &Taxon::new("f__F")),
        Err(CuratomaticError::Invariant(_))
    ));
    assert_eq!(label(&tree, 5), "90.0:f__F|RED=0.700");
    assert_eq!(label(&tree, 3), "G2|RED=1.000");
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
    let mut found = tree.tree.ids().filter_map(|node| tree.tree.label(node)).filter(|l| l.contains(needle));
    let label = found.next().unwrap();
    assert!(found.next().is_none(), "more than one label contains {needle}");
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
    tree.check_no_same_rank_taxa().unwrap();
}

#[test]
fn a_nested_taxon_is_not_moved_onto_the_node_of_its_ancestor() {
    let tree = refine_genus(NESTED_ON_ANCESTOR_NODE);
    assert_eq!(label_with(&tree, "g__B"), "50.0:g__B|RED=0.930");
    assert_eq!(label_with(&tree, "g__A"), "100.0:g__A|RED=0.880");
    tree.check_no_same_rank_taxa().unwrap();
}

#[test]
fn nested_inside_finds_the_closest_named_ancestor_of_the_same_rank() {
    let tree = tree(&nested_tree("100.0:g__A; g__")).unwrap();
    let b = tree
        .tree
        .ids()
        .find(|&node| tree.tree.label(node).is_some_and(|l| l.contains("g__B")))
        .unwrap();
    let (genus, family) = (Taxon::new("g__B"), Taxon::new("f__F"));
    assert_eq!(tree.nested_inside(b, Rank::Genus, &genus), Some(&Taxon::new("g__A")));
    assert_eq!(tree.nested_inside(b, Rank::Family, &family), None);
}

#[test]
fn a_node_with_two_named_taxa_of_the_same_rank_is_an_error() {
    let tree = tree(&nested_tree("100.0:g__A; g__C")).unwrap();
    assert!(matches!(tree.check_no_same_rank_taxa(), Err(CuratomaticError::Invariant(_))));
}

#[test]
fn a_placeholder_does_not_count_as_a_second_taxon_of_the_rank() {
    let tree = tree(&nested_tree("100.0:g__A; g__")).unwrap();
    tree.check_no_same_rank_taxa().unwrap();
}
