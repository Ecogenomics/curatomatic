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
    refine(newick, Rank::Genus)
}

/// Refine one rank of a tree, as is done for each rank when curatomatic runs.
fn refine(newick: &str, rank: Rank) -> Tree {
    let mut tree = tree(newick).unwrap();
    let min_bs = BootstrapValue::new(95.0);
    let refinement = tree
        .refine_existing_labels(rank, min_bs, &red_dict(), &Metadata::default())
        .unwrap();
    tree.process_refinement(&refinement, rank, min_bs, &red_dict())
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

/// `g__B` (bootstrap `{bs}`) sits on the genus median (RED 0.90). Its parent (RED 0.85, bootstrap
/// 100) is within the genus RED window, but further from the median.
fn further_from_median_tree(bs: &str) -> String {
    format!(
        "((((G1|RED=1.000:0.1)'s__B one|RED=1.000':0.1,(G2|RED=1.000:0.1)'s__B two|RED=1.000':0.1)\
         '{bs}:g__B|RED=0.900':0.1,(G3|RED=1.000:0.1)'s__C three|RED=1.000':0.1)'100.0|RED=0.850':0.1,\
         (G4|RED=1.000:0.1)'s__D four|RED=1.000':0.1)'100.0:f__F|RED=0.700';"
    )
}

#[test]
fn an_unsupported_taxon_is_moved_up_to_a_supported_node_further_from_the_median() {
    // Support is traded for RED: the parent is further from the median, but is supported.
    let tree = refine_genus(&further_from_median_tree("50.0"));
    assert_eq!(label_with(&tree, "g__B"), "100.0:g__B|RED=0.850");
    assert_eq!(label_with(&tree, "RED=0.900"), "50.0|RED=0.900");
}

#[test]
fn a_supported_taxon_is_not_moved_up_further_from_the_median() {
    let tree = refine_genus(&further_from_median_tree("100.0"));
    assert_eq!(label_with(&tree, "g__B"), "100.0:g__B|RED=0.900");
    assert_eq!(label_with(&tree, "RED=0.850"), "100.0|RED=0.850");
}

// ---- A taxon is never moved up above a taxon of a higher rank ----

/// The r226 bac120 `f__Sporomusaceae` case: `f__F` (bootstrap 54, RED 0.789) is on the node of
/// its order, `{taxa}`. Its parent (bootstrap 71, RED 0.749) is better supported and closer to the
/// family median, so it would otherwise qualify.
fn family_with_its_order(taxa: &str) -> String {
    format!(
        "((((G1|RED=1.000:0.1)'s__A one|RED=1.000':0.1,(G2|RED=1.000:0.1)'s__A two|RED=1.000':0.1)\
         '54.0:{taxa}|RED=0.789':0.1,(G3|RED=1.000:0.1)'s__C three|RED=1.000':0.1)'71.0|RED=0.749':0.1,\
         (G4|RED=1.000:0.1)'s__D four|RED=1.000':0.1)'100.0:c__C|RED=0.400';"
    )
}

#[test]
fn a_family_without_its_order_on_its_node_is_moved_up() {
    // Without the order the fixture qualifies, so the test below shows the order is what stops it.
    let tree = refine(&family_with_its_order("f__F"), Rank::Family);
    assert_eq!(label_with(&tree, "f__F"), "71.0:f__F|RED=0.749");
}

#[test]
fn a_family_on_the_node_of_its_order_is_not_moved_up_above_it() {
    let tree = refine(&family_with_its_order("o__O; f__F"), Rank::Family);
    assert_eq!(label_with(&tree, "f__F"), "54.0:o__O; f__F|RED=0.789");
    assert_eq!(label_with(&tree, "RED=0.749"), "71.0|RED=0.749");
}

#[test]
fn a_genus_is_not_moved_up_to_a_node_whose_other_child_holds_a_family() {
    // `further_from_median_tree("50.0")` moves g__B up to RED 0.85, unless G3 is in a family.
    let newick = further_from_median_tree("50.0").replace("'s__C three", "'f__X; s__C three");
    let tree = refine_genus(&newick);
    assert_eq!(label_with(&tree, "g__B"), "50.0:g__B|RED=0.900");
    assert_eq!(label_with(&tree, "RED=0.850"), "100.0|RED=0.850");
}

#[test]
fn a_higher_rank_on_the_candidate_node_itself_does_not_stop_the_move() {
    let newick =
        further_from_median_tree("50.0").replace("'100.0|RED=0.850'", "'100.0:f__X|RED=0.850'");
    let tree = refine_genus(&newick);
    assert_eq!(label_with(&tree, "g__B"), "100.0:f__X; g__B|RED=0.850");
}

#[test]
fn the_walk_stops_at_the_first_node_with_a_higher_rank_below() {
    // G3 is in f__X, below Y (the parent of g__B, RED 0.88), the first candidate. X above it holds
    // a family itself, which alone would not stop a move onto X, but f__X is below X too, so
    // g__B is moved to neither.
    let newick = nested_tree("100.0:f__Y").replace("'s__A three", "'f__X; s__A three");
    let tree = refine_genus(&newick);
    assert_eq!(label_with(&tree, "g__B"), "50.0:g__B|RED=0.930");
    assert_eq!(label_with(&tree, "RED=0.880"), "100.0|RED=0.880");
}

#[test]
fn nested_inside_finds_the_closest_named_ancestor_of_the_same_rank() {
    let tree = tree(&nested_tree("100.0:g__A")).unwrap();
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
//
// The file restores the labels of the dummy nodes of an old style tree: a leaf gets every named
// taxon of its row that is not on its path (see `Tree::assign_leaf_taxonomy`).

fn taxonomy_file(contents: &str) -> TaxonomyFile {
    TaxonomyFile::from_reader(contents.as_bytes(), Path::new("taxonomy.tsv")).unwrap()
}

/// A tree in the r237 format: leaves without a RED, and `|red=`. Post-order: G1 (0), G2 (1),
/// `g__A` (2), G3 (3), `f__F` (4), G4 (5), `c__C` (6), root (7).
const NEW_STYLE: &str = "(((G1:0.1,G2:0.1)'100:g__A|red=0.95':0.1,G3:0.1)'100:f__F|red=0.70':0.1,\
     (G4:0.1)'100:c__C|red=0.40':0.1)'d__D|red=0.0000';";

/// The taxonomy file for [`NEW_STYLE`]: the taxa on each genome's path, then its own. G3's genus
/// is a placeholder, and G9 is not on the tree.
const NEW_STYLE_TAXONOMY: &str = "\
G1\td__D;p__;c__;o__;f__F;g__A;s__A one
G2\td__D;p__;c__;o__;f__F;g__A;s__{unresolved} sp2
G3\td__D;p__;c__;o__;f__F;g__;s__B three
G4\td__D;p__;c__C;o__O2;f__F2;g__G;s__G four
G9\td__D;p__P;c__C;o__O;f__F;g__A;s__A nine
";

fn assigned(newick: &str, taxonomy: &str) -> (Tree, LeafAssignment) {
    let mut tree = tree(newick).unwrap();
    let assignment = tree.assign_leaf_taxonomy(&taxonomy_file(taxonomy)).unwrap();
    (tree, assignment)
}

/// The error of assigning `taxonomy` to `newick`.
fn conflict(newick: &str, taxonomy: &str) -> CuratomaticError {
    let mut tree = tree(newick).unwrap();
    tree.assign_leaf_taxonomy(&taxonomy_file(taxonomy))
        .unwrap_err()
}

/// The taxa given to the leaf of each genome, in tree order.
fn leaf_taxa_by_gid(tree: &Tree) -> Vec<(String, String)> {
    tree.newick
        .leaves()
        .map(|leaf| {
            let leaf = &tree.leaves[&leaf];
            (leaf.gid.to_string(), leaf.taxonomy.to_string())
        })
        .collect()
}

#[test]
fn a_leaf_gets_every_taxon_of_its_row_that_is_not_on_its_path() {
    let (tree, assignment) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    assert_eq!(leaf_taxa(&tree, 0), "s__A one");
    // G4 is below a class only: the order down to the species are its own.
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
}

/// A single genome genus nested in a genus on the tree (r237 G053840885).
const NESTED_GENUS: (&str, &str) = (
    "(((G1:0.1,G2:0.1)'90:red=0.97':0.1,G3:0.1)'100:g__JAAYFO01|red=0.93':0.1,\
     G4:0.1)'d__D;p__P;c__C;o__O;f__F|red=0.0';",
    "\
G1\td__D;p__P;c__C;o__O;f__F;g__JAAYFO01;g__Neosphaerochaeta;s__Neosphaerochaeta fermentans
G2\td__D;p__P;c__C;o__O;f__F;g__JAAYFO01;s__JAAYFO01 sp2
G3\td__D;p__P;c__C;o__O;f__F;g__JAAYFO01;s__JAAYFO01 sp3
G4\td__D;p__P;c__C;o__O;f__F;g__B;s__B four
",
);

#[test]
fn a_single_genome_genus_nested_in_a_tree_genus_is_restored() {
    // r237 G053840885: the tree has `g__JAAYFO01` above five genomes, and one of them is in its
    // own genus, which was on its dummy node. A node with only a bootstrap is between them.
    let (newick, taxonomy) = NESTED_GENUS;
    let (tree, _) = assigned(newick, taxonomy);
    assert_eq!(
        leaf_taxa_by_gid(&tree)[0],
        (
            "G1".to_owned(),
            "g__Neosphaerochaeta; s__Neosphaerochaeta fermentans".to_owned()
        )
    );
    assert_eq!(leaf_taxa_by_gid(&tree)[1].1, "s__JAAYFO01 sp2");
    // The tree's genus is unchanged, and the nested genus is the genome's genus.
    assert!(tree.taxon_labels.contains(&Taxon::new("g__JAAYFO01")));
    let taxonomies = tree.leaf_taxonomies().unwrap();
    assert_eq!(taxonomies[0].1[5], Taxon::new("g__Neosphaerochaeta"));
    assert_eq!(taxonomies[1].1[5], Taxon::new("g__JAAYFO01"));
}

/// A single genome family nested in a family on the tree (r226 ar53 G018920245).
const NESTED_FAMILY: (&str, &str) = (
    "((G1:0.1,G2:0.1)'100:f__JACQIB01|red=0.75':0.1,G3:0.1)'d__D;p__P;c__C;o__O|red=0.0';",
    "\
G1\td__D;p__P;c__C;o__O;f__JACQIB01;f__JAHLMS01;g__JAHLMS01;s__JAHLMS01 sp1
G2\td__D;p__P;c__C;o__O;f__JACQIB01;g__JACQIB01;s__JACQIB01 sp2
G3\td__D;p__P;c__C;o__O;f__X;g__X;s__X three
",
);

#[test]
fn a_single_genome_family_nested_in_a_tree_family_is_restored() {
    // r226 ar53 G018920245: `f__JAHLMS01; g__JAHLMS01; s__…` was its dummy node.
    let (newick, taxonomy) = NESTED_FAMILY;
    let (tree, _) = assigned(newick, taxonomy);
    assert_eq!(
        leaf_taxa_by_gid(&tree),
        [
            (
                "G1".to_owned(),
                "f__JAHLMS01; g__JAHLMS01; s__JAHLMS01 sp1".to_owned()
            ),
            ("G2".to_owned(), "g__JACQIB01; s__JACQIB01 sp2".to_owned()),
            ("G3".to_owned(), "f__X; g__X; s__X three".to_owned()),
        ]
    );
}

/// A single genome order nested in a family on the tree (r226 bac120 G017303975).
const NESTED_ORDER: (&str, &str) = (
    "((G1:0.1,G2:0.1)'100:f__JABDFU01|red=0.75':0.1,G3:0.1)'100:d__D;o__JABDFU01|red=0.6':0.1;",
    "\
G1\td__D;p__;c__;o__JABDFU01;o__JAFLBG01;f__JABDFU01;f__JAFLBG01;g__JAFLBG01;s__JAFLBG01 sp1
G2\td__D;p__;c__;o__JABDFU01;f__JABDFU01;g__JABDFU01;s__JABDFU01 sp2
G3\td__D;p__;c__;o__JABDFU01;f__Y;g__Y;s__Y three
",
);

#[test]
fn a_single_genome_order_nested_in_a_tree_family_is_restored() {
    // r226 bac120 G017303975: the dummy node held an order, which is above the tree's nearest
    // taxon (a family) in rank, so the row has it before that family. Its place in the row does
    // not matter, only that it is not on the path.
    let (newick, taxonomy) = NESTED_ORDER;
    let (tree, _) = assigned(newick, taxonomy);
    assert_eq!(
        leaf_taxa_by_gid(&tree)[0].1,
        "o__JAFLBG01; f__JAFLBG01; g__JAFLBG01; s__JAFLBG01 sp1"
    );
}

/// A family nested in a genus on the tree, out of rank order (r226 bac120 G028712565).
const OUT_OF_RANK_ORDER: (&str, &str) = (
    "(((G1:0.1,G2:0.1)'100:f__JAQTPW01|red=0.85':0.1,G3:0.1)'100:g__JAQTPY01|red=0.81':0.1,\
     G4:0.1)'d__D;p__P;c__C;o__O;f__RBG|red=0.0';",
    "\
G1\td__D;p__P;c__C;o__O;f__RBG;f__JAQTPW01;g__JAQTPY01;g__JAQTPW01;s__JAQTPW01 sp1
G2\td__D;p__P;c__C;o__O;f__RBG;f__JAQTPW01;g__JAQTPY01;s__JAQTPY01 sp2
G3\td__D;p__P;c__C;o__O;f__RBG;g__JAQTPY01;s__JAQTPY01 sp3
G4\td__D;p__P;c__C;o__O;f__RBG;g__Z;s__Z four
",
);

#[test]
fn a_taxon_nested_out_of_rank_order_on_the_path_is_not_given_to_the_leaf() {
    // r226 bac120 G028712565: the tree has a family inside a genus. The row is in rank order, so
    // the tree's genus is after the family that is nearest the leaf, but it is on the path.
    let (newick, taxonomy) = OUT_OF_RANK_ORDER;
    let (tree, _) = assigned(newick, taxonomy);
    assert_eq!(
        leaf_taxa_by_gid(&tree)
            .into_iter()
            .map(|(_, taxa)| taxa)
            .collect::<Vec<_>>(),
        [
            "g__JAQTPW01; s__JAQTPW01 sp1",
            "s__JAQTPY01 sp2",
            "s__JAQTPY01 sp3",
            "g__Z; s__Z four"
        ]
    );
}

/// The nearest node of G1 and G2 holds two ranks, and the root two taxa of one rank.
const SEVERAL_TAXA: (&str, &str) = (
    "((G1:0.1,G2:0.1)'100:f__X;g__X|red=0.9':0.1,G3:0.1)'d__D;p__A;p__B|red=0.2':0.1;",
    "\
G1\td__D;p__A;p__B;c__;o__;f__X;g__X;s__X one
G2\td__D;p__A;p__B;c__;o__;f__X;g__X;g__Y;s__Y two
G3\td__D;p__A;p__B;c__C3;o__;f__;g__;s__{unresolved} sp3
",
);

#[test]
fn every_taxon_of_a_nearest_node_with_several_taxa_is_on_the_path() {
    // The nearest node holds two ranks, and two taxa of one rank (a shared name).
    let (newick, taxonomy) = SEVERAL_TAXA;
    let (tree, _) = assigned(newick, taxonomy);
    assert_eq!(
        leaf_taxa_by_gid(&tree)
            .into_iter()
            .map(|(_, taxa)| taxa)
            .collect::<Vec<_>>(),
        ["s__X one", "g__Y; s__Y two", "c__C3; s__{unresolved} sp3"]
    );
}

#[test]
fn the_taxa_on_the_tree_are_never_changed_by_the_taxonomy_file() {
    let before = tree(NEW_STYLE).unwrap();
    let (after, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
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
fn a_leaf_with_no_taxon_on_its_path_gets_its_whole_row() {
    let (tree, _) = assigned(
        "((G1:0.1,G2:0.1)'100:red=0.5':0.1)'100|red=0.0';",
        "G1\td__D;p__P;c__C;o__O;f__F;g__G;s__G one\nG2\td__D2;s__H two\n",
    );
    assert_eq!(
        leaf_taxa(&tree, 0),
        "d__D; p__P; c__C; o__O; f__F; g__G; s__G one"
    );
}

#[test]
fn an_old_style_tree_gets_nothing_from_a_matching_taxonomy_file() {
    // The dummy nodes are still on the tree, so every taxon of each row is already on its path.
    let taxonomy = "G1\td__;f__F;g__A;s__A one\n\
                    G2\td__;f__F;s__B b\n\
                    G3\td__;f__F;g__A;s__A one\n";
    let (tree, assignment) = assigned(TWO_GENOMES, taxonomy);
    assert_eq!(assignment, LeafAssignment::default());
    assert_eq!(leaf_taxa(&tree, 0), "");
    assert_eq!(leaf_taxa(&tree, 3), "");
}

#[test]
fn a_leaf_missing_from_the_taxonomy_file_is_an_error() {
    let mut tree = tree(NEW_STYLE).unwrap();
    let file = taxonomy_file("G1\td__D;f__F;g__A;s__A one\n");
    assert_eq!(
        tree.assign_leaf_taxonomy(&file),
        Err(CuratomaticError::MissingFromTaxonomyFile {
            n: 1,
            first: GenomeId::new("G2")
        })
    );
}

#[test]
fn a_row_without_the_nearest_taxon_on_the_tree_is_an_error() {
    // G1 is in `g__A` on the tree, but `g__Other` in the file.
    let taxonomy = NEW_STYLE_TAXONOMY.replace("f__F;g__A;s__A one", "f__F;g__Other;s__A one");
    assert_eq!(
        conflict(NEW_STYLE, &taxonomy),
        CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G1"),
            message: "the tree places it in g__A, which is not in its row of the taxonomy file"
                .to_owned()
        }
    );
}

#[test]
fn a_row_missing_one_of_the_nearest_nodes_taxa_is_an_error() {
    let newick = "((G1:0.1,G2:0.1)'100:p__A;p__B|red=0.3':0.1)'100|red=0.0';";
    let taxonomy = "G1\td__;p__A;s__X one\nG2\td__;p__A;p__B;s__X two\n";
    assert_eq!(
        conflict(newick, taxonomy),
        CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G1"),
            message: "the tree places it in p__B, which is not in its row of the taxonomy file"
                .to_owned()
        }
    );
}

#[test]
fn a_tree_taxon_above_the_anchor_that_is_not_in_the_row_is_an_error() {
    // The tree puts G1 in p__A, its row says p__B. The anchor (g__X) agrees, so only a check of
    // the whole path catches it. Before this check G1 silently got p__B on its leaf.
    let newick = "((G1:0.1,G2:0.1)'100:g__X|red=0.9':0.1,G3:0.1)'d__D;p__A|red=0.0';";
    let taxonomy = "\
G1\td__D;p__B;c__;o__;f__;g__X;s__X one
G2\td__D;p__A;c__;o__;f__;g__X;s__X two
G3\td__D;p__A;c__;o__;f__;g__;s__Y three
";
    assert_eq!(
        conflict(newick, taxonomy),
        CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G1"),
            message: "the tree places it in p__A, which is not in its row of the taxonomy file"
                .to_owned()
        }
    );
}

#[test]
fn a_row_placeholder_where_the_tree_has_a_name_is_an_error() {
    // The tree puts G3 in f__F, its row has no family.
    let newick = "(((G1:0.1,G2:0.1)'100:g__X|red=0.9':0.1,G3:0.1)'100:f__F|red=0.7':0.1,G4:0.1)\
                  'd__D|red=0.0';";
    let taxonomy = "\
G1\td__D;p__;c__;o__;f__F;g__X;s__X one
G2\td__D;p__;c__;o__;f__F;g__X;s__X two
G3\td__D;p__;c__;o__;f__;g__;s__Y three
G4\td__D;p__;c__;o__;f__;g__;s__Z four
";
    assert!(matches!(
        conflict(newick, taxonomy),
        CuratomaticError::TaxonomyFileConflict { gid, message }
            if gid.as_str() == "G3" && message.contains("f__F")
    ));
}

#[test]
fn every_tree_taxon_that_is_not_in_the_row_is_named() {
    // Both the phylum and the class on the tree are missing from G1's row.
    let newick = "((G1:0.1,G2:0.1)'100:g__X|red=0.9':0.1,G3:0.1)'d__D;p__A;c__C|red=0.0';";
    let taxonomy = "\
G1\td__D;p__B;c__;o__;f__;g__X;s__X one
G2\td__D;p__A;c__C;o__;f__;g__X;s__X two
G3\td__D;p__A;c__C;o__;f__;g__;s__Y three
";
    assert_eq!(
        conflict(newick, taxonomy),
        CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G1"),
            message: "the tree places it in p__A and c__C, which are not in its row of the \
                      taxonomy file"
                .to_owned()
        }
    );
}

#[test]
fn every_genomes_full_taxonomy_matches_its_row_after_assignment() {
    // The named taxa on each genome's path and leaf are exactly those of its row, including the
    // nested and out of rank order cases.
    for (newick, taxonomy) in [
        (NEW_STYLE, NEW_STYLE_TAXONOMY),
        NESTED_GENUS,
        NESTED_FAMILY,
        NESTED_ORDER,
        OUT_OF_RANK_ORDER,
        SEVERAL_TAXA,
    ] {
        let file = taxonomy_file(taxonomy);
        let (tree, _) = assigned(newick, taxonomy);
        for (gid, taxa) in tree.full_taxonomies().unwrap() {
            let on_tree: BTreeSet<Taxon> = taxa.into_iter().filter(Taxon::is_named).collect();
            let in_row: BTreeSet<Taxon> = file.get(&gid).unwrap().iter_ranks().cloned().collect();
            assert_eq!(on_tree, in_row, "{gid} of {newick}");
        }
    }
}

#[test]
fn a_taxon_of_a_row_that_is_elsewhere_on_the_tree_is_an_error() {
    // G4 is not below `f__F`, but its row puts it there.
    let taxonomy = NEW_STYLE_TAXONOMY.replace("o__O2;f__F2", "o__O2;f__F");
    assert_eq!(
        conflict(NEW_STYLE, &taxonomy),
        CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G4"),
            message: "f__F is for this genome only, but it is already on the tree".to_owned()
        }
    );
}

#[test]
fn a_taxon_that_two_leaves_would_get_is_an_error() {
    // Both G1 and G2 have `g__N` below the tree's genus, which would need a node above both.
    let taxonomy = "\
G1\td__D;f__F;g__A;g__N;s__N one
G2\td__D;f__F;g__A;g__N;s__N two
G3\td__D;f__F;s__B three
G4\td__D;c__C;s__C four
";
    assert_eq!(
        conflict(NEW_STYLE, taxonomy),
        CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G2"),
            message: "g__N would be given to both G1 and G2".to_owned()
        }
    );
}

/// An old style tree, with a dummy node above each genome. Below `f__F` is a single genome
/// genus nested in a genus, a single genome order nested in a family, and a family nested in a
/// genus (out of rank order).
const ZOMBIE: &str = "\
((((('G1|RED=1.0':0.0,'D-G1|RED=1.0':0.0)'g__N; s__N one|RED=1.0':0.1,\
('G2|RED=1.0':0.0,'D-G2|RED=1.0':0.0)'s__A two|RED=1.0':0.1)'100.0:g__A|RED=0.95':0.1,\
('G3|RED=1.0':0.0,'D-G3|RED=1.0':0.0)'o__Q; f__Q; g__Q; s__Q three|RED=1.0':0.1)'100.0:f__F|RED=0.7':0.1,\
((('G4|RED=1.0':0.0,'D-G4|RED=1.0':0.0)'g__M; s__M four|RED=1.0':0.1,\
('G5|RED=1.0':0.0,'D-G5|RED=1.0':0.0)'s__{unresolved} sp5|RED=1.0':0.1)'99.0:f__In|RED=0.85':0.1,\
('G6|RED=1.0':0.0,'D-G6|RED=1.0':0.0)'s__Out six|RED=1.0':0.1)'100.0:g__Out|RED=0.8':0.1)\
'100.0:o__O|RED=0.5':0.1)'d__D; p__P; c__C|RED=0.0';";

/// [`ZOMBIE`] with each dummy pair collapsed into its genome.
const ZOMBIE_COLLAPSED: &str = "\
((((G1:0.1,G2:0.1)'100.0:g__A|red=0.95':0.1,G3:0.1)'100.0:f__F|red=0.7':0.1,\
((G4:0.1,G5:0.1)'99.0:f__In|red=0.85':0.1,G6:0.1)'100.0:g__Out|red=0.8':0.1)\
'100.0:o__O|red=0.5':0.1)'d__D;p__P;c__C|red=0.0';";

/// The taxa on the path of each genome of `tree` (its dummy node included), in rank order, as a
/// taxonomy file.
fn rows_from_zombie(tree: &NewickTree) -> String {
    let mut out = String::new();
    for leaf in tree.leaves() {
        let gid = NodeLabel::from_node(tree, leaf).unwrap().gid.unwrap();
        if gid.as_str().starts_with("D-") {
            continue;
        }
        let mut taxa: Vec<Taxon> = Vec::new();
        let mut current = tree.parent(leaf);
        while let Some(node) = current {
            if let Some(taxonomy) = NodeLabel::from_node(tree, node).unwrap().taxonomy {
                // Each node's taxa go in front: the root's end up first, in path order.
                let mut node_taxa: Vec<Taxon> = taxonomy.iter_ranks().cloned().collect();
                node_taxa.append(&mut taxa);
                taxa = node_taxa;
            }
            current = tree.parent(node);
        }
        // Rank order, keeping path order within a rank (a stable sort), with the placeholder of
        // each rank that has no taxon, as in the taxonomy file.
        for rank in Rank::iter() {
            if !taxa.iter().any(|taxon| taxon.rank().unwrap() == rank) {
                taxa.push(Taxon::new(format!("{}__", rank.prefix())));
            }
        }
        taxa.sort_by_key(|taxon| taxon.rank().unwrap());
        let taxa: Vec<&str> = taxa.iter().map(Taxon::as_str).collect();
        out.push_str(gid.as_str());
        out.push('\t');
        out.push_str(&taxa.join(";"));
        out.push('\n');
    }
    out
}

#[test]
fn removing_the_dummy_nodes_and_assigning_the_file_restores_their_taxa() {
    let zombie = NewickTree::parse(ZOMBIE).unwrap();
    // The taxa on each genome's dummy node.
    let dummy_taxa: Vec<(String, String)> = zombie
        .leaves()
        .filter_map(|leaf| {
            let gid = NodeLabel::from_node(&zombie, leaf).unwrap().gid.unwrap();
            let parent = zombie.parent(leaf).unwrap();
            let taxa = NodeLabel::from_node(&zombie, parent).unwrap().taxonomy;
            (!gid.as_str().starts_with("D-")).then(|| (gid.to_string(), taxa.unwrap().to_string()))
        })
        .collect();
    assert_eq!(dummy_taxa.len(), 6);

    let rows = rows_from_zombie(&zombie);
    assert!(rows.contains("G3\td__D;p__P;c__C;o__O;o__Q;f__F;f__Q;g__Q;s__Q three\n"));
    assert!(rows.contains("G6\td__D;p__P;c__C;o__O;f__;g__Out;s__Out six\n"));
    let (tree, assignment) = assigned(ZOMBIE_COLLAPSED, &rows);
    assert_eq!(leaf_taxa_by_gid(&tree), dummy_taxa);
    assert_eq!(assignment.n_leaves, 6);

    // Every taxon of the old tree is on the new one, on one node each.
    let old: BTreeSet<Taxon> = Tree::from_newick(zombie)
        .unwrap()
        .taxon_labels
        .into_iter()
        .collect();
    assert_eq!(tree.taxon_labels, old);
}

// ---- The loaded inputs, with the dummy nodes put back (`--debug`) ----

/// [`ZOMBIE_COLLAPSED`] loaded with the taxonomy file read off [`ZOMBIE`].
fn zombie_loaded() -> (NewickTree, Tree) {
    let zombie = NewickTree::parse(ZOMBIE).unwrap();
    let rows = rows_from_zombie(&zombie);
    let (tree, _) = assigned(ZOMBIE_COLLAPSED, &rows);
    (zombie, tree)
}

#[test]
fn the_zombie_tree_of_the_loaded_inputs_is_the_old_tree() {
    let (zombie, tree) = zombie_loaded();
    let rebuilt = NewickTree::parse(&tree.to_zombie_newick()).unwrap();
    // Node for node: the same shape, branch lengths, and label values. The label text differs
    // (`; ` and `;`, `RED=` and `red=`), its values never do.
    assert_eq!(rebuilt.ids().count(), zombie.ids().count());
    for (a, b) in zombie.ids().zip(rebuilt.ids()) {
        assert_eq!(
            zombie.children(a).len(),
            rebuilt.children(b).len(),
            "node {a}"
        );
        assert_eq!(zombie.length(a), rebuilt.length(b), "node {a}");
        assert_eq!(
            NodeLabel::from_node(&zombie, a).unwrap(),
            NodeLabel::from_node(&rebuilt, b).unwrap(),
            "node {a}: {:?} vs {:?}",
            zombie.label(a),
            rebuilt.label(b)
        );
    }
}

#[test]
fn the_zombie_tree_puts_the_leafs_branch_length_above_its_dummy_node() {
    let (_, tree) = zombie_loaded();
    let text = tree.to_zombie_newick();
    assert!(
        text.starts_with(
            "(((((G1|RED=1.0:0.0,D-G1|RED=1.0:0.0)'g__N; s__N one|RED=1.0':0.1,\
             (G2|RED=1.0:0.0,D-G2|RED=1.0:0.0)'s__A two|RED=1.0':0.1)'100.0:g__A|red=0.95':0.1,"
        ),
        "{text}"
    );
    // Only formatting: the tree itself is unchanged, and its leaves are still labelled G1, ...
    assert_eq!(label(&tree, 0), "G1");
}

#[test]
fn the_zombie_tree_gives_a_leaf_without_taxa_an_empty_dummy_node() {
    let tree = tree("((G1:0.1,G2:0.2)'100:g__A|red=0.9':0.3)'d__D|red=0.0';").unwrap();
    assert_eq!(
        tree.to_zombie_newick(),
        "(((G1|RED=1.0:0.0,D-G1|RED=1.0:0.0)|RED=1.0:0.1,(G2|RED=1.0:0.0,D-G2|RED=1.0:0.0)\
         |RED=1.0:0.2)'100:g__A|red=0.9':0.3)'d__D|red=0.0';\n"
    );
}

#[test]
fn full_taxonomies_are_the_taxonomy_file_that_was_loaded() {
    // Every nested taxon is kept, in rank order and path order within a rank, with placeholders.
    let (zombie, tree) = zombie_loaded();
    let mut text = String::new();
    for (gid, taxa) in tree.full_taxonomies().unwrap() {
        let taxa: Vec<&str> = taxa.iter().map(Taxon::as_str).collect();
        text.push_str(gid.as_str());
        text.push('\t');
        text.push_str(&taxa.join(";"));
        text.push('\n');
    }
    assert_eq!(text, rows_from_zombie(&zombie));
    assert!(
        text.contains("G4\td__D;p__P;c__C;o__O;f__In;g__Out;g__M;s__M four\n"),
        "{text}"
    );
}

#[test]
fn full_taxonomies_fill_a_missing_rank_with_its_placeholder() {
    let (tree, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    let rows = tree.full_taxonomies().unwrap();
    let g3: Vec<&str> = rows[2].1.iter().map(Taxon::as_str).collect();
    assert_eq!(
        g3,
        ["d__D", "p__", "c__", "o__", "f__F", "g__", "s__B three"]
    );
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
         'd__D|red=0.0';",
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
        "G1\td__D;p__P;c__C;o__O;f__F;s__{unresolved} G1\n\
         G2\td__D;p__P;c__C;o__O;f__F;s__{unresolved} G2\n\
         G3\td__D;p__P;c__C;o__O;f__F;s__{unresolved} G3\n\
         G4\td__D;p__P;c__C;o__O;f__X;g__X;s__X four\n",
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
fn a_missing_genus_above_a_named_species_is_created() {
    // Species never block a novel genus: every genome has a named species.
    let (tree, novel) = fill_novel(
        MISSING_GENUS,
        &[(0, "s__{unresolved} G1"), (1, "s__Chlorobium limicola")],
    );
    assert_eq!(novel, BTreeSet::from([Taxon::new("g__G1")]));
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__P;c__C;o__O;f__F;g__G1;s__{unresolved} G1",
            "G2\td__D;p__P;c__C;o__O;f__F;g__G1;s__Chlorobium limicola"
        ]
    );
}

#[test]
fn a_missing_family_above_a_named_genus_is_not_created() {
    // The highest family candidate (RED 0.75) is above `g__G`.
    let (tree, novel) = fill_novel(
        "(((G1:0.1,G2:0.1)'100:g__G|red=0.95':0.1)'100|red=0.75':0.1)'d__D;p__P;c__C;o__O|red=0.0';",
        &[(0, "s__G one"), (1, "s__G two")],
    );
    assert!(novel.is_empty());
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__P;c__C;o__O;f__;g__G;s__G one",
            "G2\td__D;p__P;c__C;o__O;f__;g__G;s__G two"
        ]
    );
}

#[test]
fn a_wgs_style_name_below_blocks_creation() {
    // The highest family candidate (RED 0.75) is above `g__CAYYKD01`, so no family is created.
    // G3's genus is decided on its own: its only candidate is the leaf itself.
    let (_, novel) = fill_novel(
        "(((G1:0.1,G2:0.1)'100:g__CAYYKD01|red=0.95':0.1,G3:0.1)'100|red=0.75':0.1)\
         'd__D;p__P;c__C;o__O|red=0.0';",
        &[
            (0, "s__{unresolved} G1"),
            (1, "s__{unresolved} G2"),
            (3, "s__{unresolved} G3"),
        ],
    );
    assert_eq!(novel, BTreeSet::from([Taxon::new("g__G3")]));
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
fn the_lower_ranks_of_a_skipped_genome_are_still_filled() {
    // The phylum is skipped at X (RED 0.2), above `c__C`. Nothing named is below Z (RED 0.75), the
    // highest family candidate, so the family is filled. Z is too far from the genus median to be
    // a genus candidate, so each genome gets its own genus.
    let (tree, novel) = fill_novel(
        "((((G1:0.1,G2:0.1)'100|red=0.75':0.1)'100:c__C;o__O|red=0.45':0.1)'100|red=0.2':0.1)\
         'd__D|red=0.0';",
        &[(0, "s__{unresolved} G1"), (1, "s__{unresolved} G2")],
    );
    assert_eq!(
        novel,
        BTreeSet::from([
            Taxon::new("f__G1"),
            Taxon::new("g__G1"),
            Taxon::new("g__G2")
        ])
    );
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__;c__C;o__O;f__G1;g__G1;s__{unresolved} G1",
            "G2\td__D;p__;c__C;o__O;f__G1;g__G2;s__{unresolved} G2"
        ]
    );
}

#[test]
fn a_lower_rank_of_a_skipped_genome_is_skipped_on_its_own_candidate() {
    // The phylum is skipped above `c__C`, and the family is skipped too, because its own highest
    // candidate Z (RED 0.75) is above `g__G`. Only the phylum and family are left empty.
    let (tree, novel) = fill_novel(
        "((((G1:0.1,G2:0.1)'100:g__G|red=0.95':0.1)'100|red=0.75':0.1)'100:c__C;o__O|red=0.45':0.1)\
         'd__D|red=0.0';",
        &[(0, "s__G one"), (1, "s__G two")],
    );
    assert!(novel.is_empty());
    assert_eq!(
        written_taxonomies(&tree),
        [
            "G1\td__D;p__;c__C;o__O;f__;g__G;s__G one",
            "G2\td__D;p__;c__C;o__O;f__;g__G;s__G two"
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

#[test]
fn taxon_counts_include_the_taxa_of_internal_nodes_and_leaves() {
    let tree = tree_with_leaf_taxa(
        "((G1:0.1,G2:0.1)'100:f__F; g__A|red=0.9':0.1,G3:0.1)'d__D|red=0.0';",
        &[(0, "s__A one"), (1, "s__A two"), (3, "g__B; s__B three")],
    )
    .unwrap();
    let counts = tree.taxon_counts().unwrap();
    assert_eq!(
        counts
            .iter()
            .map(|(taxon, n)| (taxon.as_str(), n))
            .collect::<Vec<_>>(),
        [
            ("d__D", 1),
            ("f__F", 1),
            ("g__A", 1),
            ("g__B", 1),
            ("s__A one", 1),
            ("s__A two", 1),
            ("s__B three", 1),
        ]
    );
}

// ---- Every leaf holds its species ----

#[test]
fn every_leaf_holds_its_species_once_the_taxonomy_file_is_assigned() {
    let (tree, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    assert_eq!(tree.check_leaves_have_species(), Ok(()));
}

#[test]
fn a_leaf_without_any_taxon_is_an_error() {
    let tree = tree(NEW_STYLE).unwrap();
    assert_eq!(
        tree.check_leaves_have_species(),
        Err(CuratomaticError::LeafWithoutSpecies {
            n: 4,
            first: GenomeId::new("G1"),
        })
    );
}

#[test]
fn a_leaf_with_taxa_but_no_species_is_an_error() {
    // G4 holds an order, family and genus, but its species is a placeholder.
    let taxonomy = NEW_STYLE_TAXONOMY.replace("s__G four", "s__");
    let (tree, _) = assigned(NEW_STYLE, &taxonomy);
    assert_eq!(leaf_taxa(&tree, 5), "o__O2; f__F2; g__G");
    assert_eq!(
        tree.check_leaves_have_species(),
        Err(CuratomaticError::LeafWithoutSpecies {
            n: 1,
            first: GenomeId::new("G4"),
        })
    );
}

#[test]
fn a_species_on_an_internal_node_is_not_the_leafs() {
    // The species is on the node above G4, so G4 itself does not hold it.
    let newick = NEW_STYLE.replace("'100:c__C|red=0.40'", "'100:c__C;s__G four|red=0.40'");
    let (tree, _) = assigned(&newick, NEW_STYLE_TAXONOMY);
    assert_eq!(
        tree.check_leaves_have_species(),
        Err(CuratomaticError::LeafWithoutSpecies {
            n: 1,
            first: GenomeId::new("G4"),
        })
    );
}

// ---- Every leaf has a named domain ----

#[test]
fn every_leaf_has_a_domain_once_the_taxonomy_file_is_assigned() {
    let (tree, _) = assigned(NEW_STYLE, NEW_STYLE_TAXONOMY);
    assert_eq!(tree.check_leaves_have_domain(), Ok(()));
}

#[test]
fn a_leaf_without_a_named_domain_is_an_error() {
    // No domain on the tree, and G1 and G2 have a placeholder in their rows. G3 holds its own.
    let newick = "((G1:0.1,G2:0.1)'100:g__A|red=0.9':0.1,G3:0.1)'100:p__P|red=0.0';";
    let taxonomy = "\
G1\td__;p__P;c__;o__;f__;g__A;s__A one
G2\td__;p__P;c__;o__;f__;g__A;s__A two
G3\td__D;p__P;c__;o__;f__;g__B;s__B three
";
    let (tree, _) = assigned(newick, taxonomy);
    assert_eq!(leaf_taxa(&tree, 3), "d__D; g__B; s__B three");
    assert_eq!(
        tree.check_leaves_have_domain(),
        Err(CuratomaticError::LeafWithoutDomain {
            n: 2,
            first: GenomeId::new("G1"),
        })
    );
}

#[test]
fn a_domain_on_an_internal_node_counts_for_the_leaves_below_it() {
    let newick = "(((G1:0.1,G2:0.1)'100:d__D;g__A|red=0.9':0.1,G3:0.1)'100:p__P|red=0.2':0.1)\
                  '100|red=0.0';";
    let taxonomy = "\
G1\td__D;p__P;c__;o__;f__;g__A;s__A one
G2\td__D;p__P;c__;o__;f__;g__A;s__A two
G3\td__;p__P;c__;o__;f__;g__B;s__B three
";
    let (tree, _) = assigned(newick, taxonomy);
    assert_eq!(
        tree.check_leaves_have_domain(),
        Err(CuratomaticError::LeafWithoutDomain {
            n: 1,
            first: GenomeId::new("G3"),
        })
    );
}

/// `g__Alpha` is on an unsupported node (bootstrap 91, RED 0.940, inside the genus window of
/// [`red_dict`]) above a supported node (bootstrap 100, RED 0.970) with G1 and G2, and the leaf
/// G3. G4 is outside the genus. Post-order: G1 (0), G2 (1), their node (2), G3 (3), the node of
/// `g__Alpha` (4), G4 (5), root (6).
///
/// G1, G2 and G4 are in the prior taxonomy; G3 is in it only if `g3_prior` is given. The species
/// of G3 is `g3_species`.
fn genus_above_a_new_genome(g3_species: &str, g3_prior: Option<&str>) -> (Tree, Vec<Move>) {
    genus_above_a_new_genome_at("g__Alpha", "0.940", g3_species, g3_prior)
}

/// [`genus_above_a_new_genome`] with the genus named `genus` on a node with RED `red`.
fn genus_above_a_new_genome_at(
    genus: &str,
    red: &str,
    g3_species: &str,
    g3_prior: Option<&str>,
) -> (Tree, Vec<Move>) {
    let mut tree = tree_with_leaf_taxa(
        &format!(
            "(((G1:0.1,G2:0.1)'100:red=0.97':0.1,G3:0.1)'91:{genus}|red={red}':0.1,G4:0.1)\
             '100:f__F|red=0.70';"
        ),
        &[
            (0, "s__Alpha one"),
            (1, "s__Alpha two"),
            (3, g3_species),
            (5, "g__C;s__C four"),
        ],
    )
    .unwrap();
    let mut metadata = Metadata::default();
    let prior = [
        ("G1", Some(format!("f__F;{genus};s__Alpha one"))),
        ("G2", Some(format!("f__F;{genus};s__Alpha two"))),
        ("G3", g3_prior.map(str::to_owned)),
        ("G4", Some("f__F;g__C;s__C four".to_owned())),
    ];
    for (gid, taxonomy) in prior {
        if let Some(taxonomy) = taxonomy {
            metadata.set_taxonomy(
                GenomeId::new(gid),
                Taxonomy::from_string(&taxonomy).unwrap(),
            );
        }
    }
    let min_bs = BootstrapValue::new(95.0);
    let refinement = tree
        .refine_existing_labels(Rank::Genus, min_bs, &red_dict(), &metadata)
        .unwrap();
    let moves = tree
        .process_refinement(&refinement, Rank::Genus, min_bs, &red_dict())
        .unwrap();
    tree.reload().unwrap();
    (tree, moves)
}

#[test]
fn a_genus_is_not_moved_down_away_from_a_new_genome_that_its_species_names() {
    // G3 is new (not in the prior), and its species puts it in `g__Alpha`: the genus is not split.
    let (tree, moves) = genus_above_a_new_genome("s__Alpha sp3", None);
    assert_eq!(moves, []);
    assert_eq!(label(&tree, 4), "91:g__Alpha|red=0.940");
    assert_eq!(label(&tree, 2), "100:red=0.97");
}

#[test]
fn a_genus_is_moved_down_away_from_a_new_genome_whose_species_names_another_genus() {
    for species in ["s__B sp3", "s__{unresolved} Alpha sp3"] {
        let (tree, moves) = genus_above_a_new_genome(species, None);
        assert_eq!(moves.len(), 1, "{species}");
        assert_eq!(label(&tree, 2), "100.0:g__Alpha|RED=0.970", "{species}");
        assert_eq!(label(&tree, 4), "91.0|RED=0.940", "{species}");
    }
}

#[test]
fn the_species_of_a_genome_in_the_prior_taxonomy_does_not_override_its_prior_genus() {
    // G3 was in the prior release, in `g__B`: its species is not used.
    let (tree, moves) = genus_above_a_new_genome("s__Alpha sp3", Some("f__F;g__B;s__Alpha sp3"));
    assert_eq!(moves.len(), 1);
    assert_eq!(label(&tree, 2), "100.0:g__Alpha|RED=0.970");
}

#[test]
fn a_genus_with_a_placeholder_name_is_moved_down_away_from_a_new_genome_its_species_names() {
    // Only a Latin genus name keeps its new genomes: `s__JANWYU01 sp3` echoes the decoration.
    let (tree, moves) =
        genus_above_a_new_genome_at("g__JANWYU01", "0.940", "s__JANWYU01 sp3", None);
    assert_eq!(moves.len(), 1);
    assert_eq!(label(&tree, 2), "100.0:g__JANWYU01|RED=0.970");
}

#[test]
fn a_genus_outside_the_red_window_is_moved_down_away_from_a_new_genome_its_species_names() {
    // The genus window of `red_dict` is 0.80 to 0.95: the species of a new genome only keeps a
    // genus whose node is inside it.
    for red in ["0.958", "0.790"] {
        let (tree, moves) = genus_above_a_new_genome_at("g__Alpha", red, "s__Alpha sp3", None);
        assert_eq!(moves.len(), 1, "{red}");
        assert_eq!(label(&tree, 2), "100.0:g__Alpha|RED=0.970", "{red}");
    }
}
