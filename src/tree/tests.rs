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
