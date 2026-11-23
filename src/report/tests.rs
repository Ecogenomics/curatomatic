use std::io;

use super::*;
use crate::genome::gid_to_canonical;
use crate::taxonomy::Taxonomy;
use crate::tree::fill_in_missing_ranks;

/// Genus: parent threshold 0.8, median 0.9, child threshold 0.95, so `q < 0.5` is RED < 0.85.
fn red_dict() -> Red {
    Red::from_tsv("rank\tmedian\nPhylum\t0.3\nClass\t0.4\nOrder\t0.5\nFamily\t0.7\nGenus\t0.9\n")
        .unwrap()
}

/// Metadata from `accession`, `ncbi_wgs_formatted` and `gtdb_taxonomy` columns, where the taxonomy
/// is the prior taxonomy (`none` if the genome has none).
fn metadata(contents: &str) -> Metadata {
    let path = Path::new("names.tsv");
    let mut meta = Metadata::from_reader(
        io::Cursor::new(contents),
        path,
        "accession",
        "ncbi_wgs_formatted",
    )
    .unwrap();
    for line in contents.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols[2] != "none" {
            let gid = gid_to_canonical(cols[0]);
            meta.set_taxonomy(gid, Taxonomy::from_string(cols[2]).unwrap());
        }
    }
    meta
}

const HEAD: &str = "accession\tncbi_wgs_formatted\tgtdb_taxonomy\n";
const ABOVE_GENUS: &str = "d__D;p__P;c__C;o__O;f__F";
const ROOT: &str = "'100.0:d__D; p__P; c__C; o__O; f__F|RED=0.000'";

/// Run curatomatic on a tree as `main` does, and return the report.
fn report(newick: &str, metadata: &Metadata) -> Vec<ReportRow> {
    let min_bs = BootstrapValue::new(95.0);
    let red_dict = red_dict();
    let mut tree = Tree::from_newick(NewickTree::parse(newick).unwrap()).unwrap();

    let mut moves = Vec::new();
    for rank in [
        Rank::Phylum,
        Rank::Class,
        Rank::Order,
        Rank::Family,
        Rank::Genus,
    ] {
        let refinement = tree
            .refine_existing_labels(rank, min_bs, &red_dict, metadata)
            .unwrap();
        moves.extend(
            tree.process_refinement(&refinement, rank, min_bs, &red_dict)
                .unwrap(),
        );
        tree.reload().unwrap();
    }

    let missing = tree.get_leaf_nodes_missing_rank().unwrap();
    let known = tree.taxon_labels.clone();
    let novel =
        fill_in_missing_ranks(&mut tree, &missing, min_bs, &red_dict, metadata, known).unwrap();
    tree.reload().unwrap();

    build_report(&ReportInput {
        tree: &tree,
        moves: &moves,
        novel: &novel,
        metadata,
        prior_taxa: &metadata.all_taxa(),
        min_bs,
        red_dict: &red_dict,
    })
    .unwrap()
}

fn row<'a>(rows: &'a [ReportRow], taxon: &str) -> &'a ReportRow {
    rows.iter()
        .find(|row| row.taxon.as_str() == taxon)
        .unwrap_or_else(|| panic!("no row for {taxon}"))
}

fn flags(rows: &[ReportRow], taxon: &str) -> Vec<Flag> {
    row(rows, taxon).flags.iter().copied().collect()
}

fn taxa(names: &[&str]) -> BTreeSet<Taxon> {
    names.iter().copied().map(Taxon::new).collect()
}

/// `g__B` is on an unsupported node over G1, below a supported node over G1 to G5.
const LARGE_GAIN_TREE: &str = "(((G1|RED=1.000)'s__B one|RED=1.000')'50.0:g__B|RED=0.930',\
     (G2|RED=1.000)'s__B two|RED=1.000',(G3|RED=1.000)'s__B three|RED=1.000',\
     (G4|RED=1.000)'s__B four|RED=1.000',(G5|RED=1.000)'s__B five|RED=1.000')\
     '100.0|RED=0.880'";

#[test]
fn a_move_up_that_gains_five_times_the_genomes_is_a_large_gain() {
    let rows = report(&format!("({LARGE_GAIN_TREE}){ROOT};"), &Metadata::default());
    let b = row(&rows, "g__B");
    assert_eq!(b.action, Action::MovedUp);
    assert_eq!(flags(&rows, "g__B"), vec![Flag::LargeGain]);
    let before = b.before.unwrap();
    assert_eq!(
        (before.n_genomes, before.bootstrap, before.red),
        (1, Some(BootstrapValue::new(50.0)), RedValue::new(0.93))
    );
    assert_eq!((b.after.n_genomes, b.after.red), (5, RedValue::new(0.88)));
    assert_eq!(rows.len(), 1, "only g__B was changed");
}

#[test]
fn a_taxon_that_was_not_changed_has_no_row_even_if_it_would_be_flagged() {
    // g__B is nested in g__A on an unsupported node, but neither moves.
    let newick = format!(
        "((((G1|RED=1.000)'s__B one|RED=1.000',(G2|RED=1.000)'s__B two|RED=1.000')\
         '50.0:g__B|RED=0.930',(G4|RED=1.000)'s__A four|RED=1.000')'100.0:g__A|RED=0.880'){ROOT};"
    );
    let rows = report(&newick, &Metadata::default());
    assert!(rows.is_empty(), "{rows:?}");
}

#[test]
fn a_taxon_moved_down_inside_another_of_the_same_rank_is_nested() {
    // The prior taxonomy puts only G1 in g__B, so g__B moves down onto its species node.
    let newick = format!(
        "((((G1|RED=1.000)'s__B one|RED=1.000',(G2|RED=1.000)'s__B two|RED=1.000')\
         '50.0:g__B|RED=0.930',(G4|RED=1.000)'s__A four|RED=1.000')'100.0:g__A|RED=0.880'){ROOT};"
    );
    let meta = metadata(&format!(
        "{HEAD}G1\tnone\t{ABOVE_GENUS};g__B;s__B one\nG2\tnone\tnone\nG4\tnone\tnone\n"
    ));
    let rows = report(&newick, &meta);
    let b = row(&rows, "g__B");
    assert_eq!(b.action, Action::MovedDown);
    assert_eq!(flags(&rows, "g__B"), vec![Flag::NestedInSameRank]);
    assert_eq!(b.nested_in, Some(Taxon::new("g__A")));
    assert!(
        rows.iter().all(|row| row.taxon.as_str() != "g__A"),
        "g__A was not changed"
    );
}

/// G1 and G2 have no genus, and sit below a supported node at RED 0.84 (q = 0.4). Their species
/// are unresolved, so a novel genus may be created above them. G3 keeps the prior taxa above genus
/// polyphyletic at the root, so they are not moved down.
fn novel_tree() -> String {
    format!(
        "(((G1|RED=1.000)'s__{{unresolved}} G1|RED=1.000',\
         (G2|RED=1.000)'s__{{unresolved}} G2|RED=1.000')'100.0|RED=0.840',\
         (G3|RED=1.000)'g__C; s__C three|RED=1.000'){ROOT};"
    )
}

#[test]
fn a_novel_taxon_of_former_members_of_a_genus_is_displaced_from_prior() {
    let meta = metadata(&format!(
        "{HEAD}G1\tABC01\t{ABOVE_GENUS};g__Old;s__Old one\nG2\tnone\tnone\n\
         G3\tnone\t{ABOVE_GENUS};g__C;s__C three\n"
    ));
    let rows = report(&novel_tree(), &meta);
    let novel = row(&rows, "g__ABC01");
    assert_eq!(novel.action, Action::Novel);
    assert_eq!(novel.before, None);
    assert_eq!(novel.after.n_genomes, 2);
    assert_eq!(
        flags(&rows, "g__ABC01"),
        vec![Flag::NovelLowRed, Flag::DisplacedFromPrior]
    );
    assert_eq!(novel.prior_taxa, taxa(&["g__Old"]));
}

#[test]
fn a_novel_taxon_of_new_genomes_is_not_displaced() {
    let meta = metadata(&format!(
        "{HEAD}G1\tABC01\tnone\nG2\tnone\tnone\nG3\tnone\t{ABOVE_GENUS};g__C;s__C three\n"
    ));
    let rows = report(&novel_tree(), &meta);
    assert_eq!(flags(&rows, "g__ABC01"), vec![Flag::NovelLowRed]);
    assert!(row(&rows, "g__ABC01").prior_taxa.is_empty());
}

/// G1 and G2 each need a novel genus (on the node above their species, RED 0.9), but the genus
/// named after each genome is already on the tree, around G3 and G4.
fn unresolvable_tree() -> String {
    format!(
        "(((G1|RED=1.000)'s__{{unresolved}} G1|RED=1.000')'100.0|RED=0.900',\
         ((G2|RED=1.000)'s__{{unresolved}} G2|RED=1.000')'100.0|RED=0.900',\
         (G3|RED=1.000)'g__G1; s__X three|RED=1.000',\
         (G4|RED=1.000)'g__G2; s__Y four|RED=1.000'){ROOT};"
    )
}

#[test]
fn every_unresolvable_novel_taxon_has_a_row_with_its_namesakes() {
    let rows = report(&unresolvable_tree(), &Metadata::default());
    assert_eq!(rows.len(), 2, "{rows:#?}");
    for (row, namesake) in rows.iter().zip(["g__G1", "g__G2"]) {
        assert_eq!(row.taxon, Taxon::new("g__UNRESOLVABLE"));
        assert_eq!(row.rank, Rank::Genus);
        assert_eq!(row.action, Action::Unresolvable);
        assert_eq!(row.flags, BTreeSet::from([Flag::NamesakeOutsideTaxon]));
        assert_eq!(row.namesakes, taxa(&[namesake]));
        assert_eq!(row.before, None);
        assert_eq!(
            (row.after.n_genomes, row.after.red),
            (1, RedValue::new(0.9))
        );
    }
}

#[test]
fn an_unresolvable_novel_taxon_without_a_row_is_an_error() {
    let tree = Tree::from_newick(NewickTree::parse(&unresolvable_tree()).unwrap()).unwrap();
    let novel = NovelTaxa {
        created: BTreeSet::new(),
        unresolvable: vec![UnresolvableTaxon {
            rank: Rank::Genus,
            node: tree.newick().ids().next().unwrap(),
            namesakes: taxa(&["g__G1"]),
        }],
    };
    let input = ReportInput {
        tree: &tree,
        moves: &[],
        novel: &novel,
        metadata: &Metadata::default(),
        prior_taxa: &BTreeSet::new(),
        min_bs: BootstrapValue::new(95.0),
        red_dict: &red_dict(),
    };
    assert_eq!(
        check_every_change_has_a_row(&[], &input),
        Err(CuratomaticError::Invariant(
            "the report has no row for these changed taxa: [g__UNRESOLVABLE (g__G1)]".to_owned()
        ))
    );
}

#[test]
fn a_novel_taxon_is_never_new_in_decoration() {
    let meta = metadata(&format!(
        "{HEAD}G1\tABC01\tnone\nG2\tnone\tnone\nG3\tnone\t{ABOVE_GENUS};g__Other;s__C three\n"
    ));
    let rows = report(&novel_tree(), &meta);
    assert!(!flags(&rows, "g__ABC01").contains(&Flag::NewInDecoration));
}

#[test]
fn a_moved_taxon_missing_from_the_prior_taxonomy_is_new_in_decoration() {
    let meta = metadata(&format!(
        "{HEAD}G1\tnone\t{ABOVE_GENUS};g__X;s__X one\nG2\tnone\tnone\nG3\tnone\tnone\n\
         G4\tnone\tnone\nG5\tnone\tnone\n"
    ));
    let rows = report(&format!("({LARGE_GAIN_TREE}){ROOT};"), &meta);
    assert_eq!(row(&rows, "g__B").action, Action::MovedUp);
    assert_eq!(
        flags(&rows, "g__B"),
        vec![Flag::LargeGain, Flag::NewInDecoration]
    );
}

/// `g__B` moves down from an unsupported node (RED `red`) onto the species node of G1, the only
/// genome the prior taxonomy puts in it. Its closest sister is `sister`, the other child of that
/// node, which joins it at RED `red`. (A genome left behind in `g__B`'s old node would get a novel
/// genus, which would be the closest sister instead.)
fn moved_down_with_sister(red: &str, sister: &str) -> Vec<ReportRow> {
    let newick =
        format!("(((G1|RED=1.000)'s__B one|RED=1.000',{sister})'50.0:g__B|RED={red}'){ROOT};");
    let meta = metadata(&format!(
        "{HEAD}G1\tnone\t{ABOVE_GENUS};g__B;s__B one\nG2\tnone\tnone\nG3\tnone\tnone\n\
         G4\tnone\tnone\n"
    ));
    let rows = report(&newick, &meta);
    assert_eq!(row(&rows, "g__B").action, Action::MovedDown);
    rows
}

#[test]
fn a_sister_that_joins_past_the_child_threshold_is_flagged() {
    let sister = "(G3|RED=1.000)'g__S; s__S three|RED=1.000'";
    let rows = moved_down_with_sister("0.960", sister);
    assert_eq!(flags(&rows, "g__B"), vec![Flag::SisterJoinPastChildThresh]);
    assert_eq!(row(&rows, "g__B").sister_taxa, taxa(&["g__S"]));
    assert!(
        rows.iter().all(|row| row.taxon.as_str() != "g__S"),
        "g__S was not changed"
    );

    let rows = moved_down_with_sister("0.940", sister);
    assert!(
        flags(&rows, "g__B").is_empty(),
        "a join below the child threshold is not flagged"
    );
}

#[test]
fn the_sisters_are_the_highest_taxa_of_the_rank_under_the_other_child() {
    // g__D is nested in g__C, so only g__C is the sister of g__B.
    let sister = "((G3|RED=1.000)'g__D; s__D three|RED=1.000',\
                  (G4|RED=1.000)'s__C four|RED=1.000')'100.0:g__C|RED=0.990'";
    let rows = moved_down_with_sister("0.960", sister);
    assert_eq!(row(&rows, "g__B").sister_taxa, taxa(&["g__C"]));
}

#[test]
fn a_change_without_a_row_is_an_error() {
    let tree = Tree::from_newick(
        NewickTree::parse(&format!("((G1|RED=1.000)'s__A one|RED=1.000'){ROOT};")).unwrap(),
    )
    .unwrap();
    let novel = NovelTaxa {
        created: taxa(&["g__Missing"]),
        unresolvable: Vec::new(),
    };
    let result = build_report(&ReportInput {
        tree: &tree,
        moves: &[],
        novel: &novel,
        metadata: &Metadata::default(),
        prior_taxa: &BTreeSet::new(),
        min_bs: BootstrapValue::new(95.0),
        red_dict: &red_dict(),
    });
    assert!(
        matches!(&result, Err(CuratomaticError::Invariant(msg)) if msg.contains("g__Missing")),
        "{result:?}"
    );
}

#[test]
fn the_tsv_has_a_header_and_leaves_unset_values_empty() {
    let rows = [ReportRow {
        taxon: Taxon::new("g__ABC01"),
        rank: Rank::Genus,
        action: Action::Novel,
        flags: [Flag::DisplacedFromPrior, Flag::NovelLowRed]
            .into_iter()
            .collect(),
        before: None,
        after: NodeStats {
            n_genomes: 2,
            bootstrap: Some(BootstrapValue::new(100.0)),
            red: RedValue::new(0.84),
        },
        prior_taxa: taxa(&["g__Old", "g__Older"]),
        nested_in: None,
        sister_taxa: BTreeSet::new(),
        namesakes: BTreeSet::new(),
    }];
    assert_eq!(
        to_tsv(&rows),
        "taxon\trank\taction\tflag\tn_genomes_before\tn_genomes_after\tbootstrap_before\t\
         bootstrap_after\tred_before\tred_after\tprior_taxon\tnested_in\tsister_taxa\t\
         namesakes\n\
         g__ABC01\tg\tnovel\tnovel_low_red,displaced_from_prior\t\t2\t\t100.0\t\t0.840\t\
         g__Old,g__Older\t\t\t\n"
    );
}

#[test]
fn the_tsv_of_an_unresolvable_taxon_lists_its_namesakes() {
    let rows = report(&unresolvable_tree(), &Metadata::default());
    let tsv = to_tsv(&rows[..1]);
    assert_eq!(
        tsv.lines().nth(1),
        Some(
            "g__UNRESOLVABLE\tg\tunresolvable\tnamesake_outside_taxon\t\t1\t\t100.0\t\t0.900\t\
             \t\t\tg__G1"
        )
    );
}
