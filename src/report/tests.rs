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
    let taxa_before = tree.taxon_nodes();

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
        taxa_before: &taxa_before,
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

#[test]
fn a_move_up_that_gains_five_times_the_genomes_is_a_large_gain() {
    let newick = format!(
        "((((G1|RED=1.000)'s__B one|RED=1.000')'50.0:g__B|RED=0.930',\
         (G2|RED=1.000)'s__B two|RED=1.000',(G3|RED=1.000)'s__B three|RED=1.000',\
         (G4|RED=1.000)'s__B four|RED=1.000',(G5|RED=1.000)'s__B five|RED=1.000')\
         '100.0|RED=0.880'){ROOT};"
    );
    let rows = report(&newick, &Metadata::default());
    let b = row(&rows, "g__B");
    assert_eq!(b.action, Some(Action::MovedUp));
    assert_eq!(flags(&rows, "g__B"), vec![Flag::LargeGain]);
    let before = b.before.unwrap();
    assert_eq!(
        (before.n_genomes, before.bootstrap, before.red),
        (1, Some(BootstrapValue::new(50.0)), RedValue::new(0.93))
    );
    assert_eq!((b.after.n_genomes, b.after.red), (5, RedValue::new(0.88)));
    assert_eq!(rows.len(), 1, "only g__B was changed or flagged");
}

#[test]
fn a_taxon_inside_another_of_the_same_rank_is_nested() {
    let newick = format!(
        "((((G1|RED=1.000)'s__B one|RED=1.000',(G2|RED=1.000)'s__B two|RED=1.000')\
         '50.0:g__B|RED=0.930',(G4|RED=1.000)'s__A four|RED=1.000')'100.0:g__A|RED=0.880'){ROOT};"
    );
    let rows = report(&newick, &Metadata::default());
    let b = row(&rows, "g__B");
    assert_eq!(b.action, None);
    assert_eq!(
        flags(&rows, "g__B"),
        vec![Flag::NestedInSameRank, Flag::UnsupportedLabel]
    );
    assert_eq!(b.nested_in, Some(Taxon::new("g__A")));
    assert!(
        rows.iter().all(|row| row.taxon.as_str() != "g__A"),
        "g__A is not flagged"
    );
}

/// G1 and G2 have no genus, and sit below a supported node at RED 0.84 (q = 0.4). G3 keeps the
/// prior taxa above genus polyphyletic at the root, so they are not moved down.
fn novel_tree() -> String {
    format!(
        "(((G1|RED=1.000)'s__X one|RED=1.000',(G2|RED=1.000)'s__X two|RED=1.000')'100.0|RED=0.840',\
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
    assert_eq!(novel.action, Some(Action::Novel));
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

#[test]
fn a_taxon_missing_from_the_prior_taxonomy_is_new_in_decoration() {
    let meta = metadata(&format!(
        "{HEAD}G1\tABC01\tnone\nG2\tnone\tnone\nG3\tnone\t{ABOVE_GENUS};g__Other;s__C three\n"
    ));
    let rows = report(&novel_tree(), &meta);
    assert_eq!(flags(&rows, "g__C"), vec![Flag::NewInDecoration]);
    assert!(
        !flags(&rows, "g__ABC01").contains(&Flag::NewInDecoration),
        "novel taxa are not"
    );
}

#[test]
fn nothing_is_new_in_decoration_without_a_prior_taxonomy() {
    let meta = metadata(&format!(
        "{HEAD}G1\tABC01\tnone\nG2\tnone\tnone\nG3\tnone\tnone\n"
    ));
    let rows = report(&novel_tree(), &meta);
    assert!(
        rows.iter()
            .all(|row| !row.flags.contains(&Flag::NewInDecoration))
    );
}

#[test]
fn sisters_that_join_past_the_child_threshold_are_flagged_together() {
    let join = |red: &str| {
        format!(
            "(((G1|RED=1.000)'g__A; s__A a|RED=1.000',(G2|RED=1.000)'g__B; s__B b|RED=1.000')\
             '100.0|RED={red}'){ROOT};"
        )
    };
    let rows = report(&join("0.970"), &Metadata::default());
    assert_eq!(flags(&rows, "g__A"), vec![Flag::SisterJoinPastChildThresh]);
    assert_eq!(row(&rows, "g__A").sister_taxa, taxa(&["g__B"]));
    assert_eq!(row(&rows, "g__B").sister_taxa, taxa(&["g__A"]));

    let rows = report(&join("0.940"), &Metadata::default());
    assert!(
        rows.is_empty(),
        "a join below the child threshold is not flagged"
    );
}

#[test]
fn the_sisters_are_the_highest_taxa_of_the_rank_under_the_other_child() {
    // g__B is nested in g__C, so only g__C is the sister of g__A.
    let newick = format!(
        "(((G1|RED=1.000)'g__A; s__A a|RED=1.000',((G2|RED=1.000)'g__B; s__B b|RED=1.000',\
         (G3|RED=1.000)'s__C c|RED=1.000')'100.0:g__C|RED=0.990')'100.0|RED=0.970'){ROOT};"
    );
    let rows = report(&newick, &Metadata::default());
    assert_eq!(row(&rows, "g__A").sister_taxa, taxa(&["g__C"]));
    assert_eq!(flags(&rows, "g__B"), vec![Flag::NestedInSameRank]);
}

/// `g__A` is on a supported node at RED 0.84 (q = 0.4) over G1 and G2. G3 keeps the prior taxa above
/// genus polyphyletic at the root.
fn split_tree() -> String {
    format!(
        "(((G1|RED=1.000)'s__A one|RED=1.000',(G2|RED=1.000)'s__A two|RED=1.000')'100.0:g__A|RED=0.840',\
         (G3|RED=1.000)'g__C; s__C three|RED=1.000'){ROOT};"
    )
}

#[test]
fn a_low_taxon_whose_prior_members_are_under_one_child_is_a_split_candidate() {
    let meta = metadata(&format!(
        "{HEAD}G1\tnone\t{ABOVE_GENUS};g__A;s__A one\nG2\tnone\tnone\n\
         G3\tnone\t{ABOVE_GENUS};g__C;s__C three\n"
    ));
    let rows = report(&split_tree(), &meta);
    assert_eq!(row(&rows, "g__A").action, None);
    assert_eq!(flags(&rows, "g__A"), vec![Flag::SplitCandidate]);
}

#[test]
fn a_taxon_whose_other_child_holds_a_prior_genome_is_not_a_split_candidate() {
    let meta = metadata(&format!(
        "{HEAD}G1\tnone\t{ABOVE_GENUS};g__A;s__A one\nG2\tnone\t{ABOVE_GENUS};g__Z;s__Z two\n\
         G3\tnone\t{ABOVE_GENUS};g__C;s__C three\n"
    ));
    let rows = report(&split_tree(), &meta);
    assert!(rows.iter().all(|row| row.taxon.as_str() != "g__A"));
}

#[test]
fn the_tsv_has_a_header_and_leaves_unset_values_empty() {
    let rows = [ReportRow {
        taxon: Taxon::new("g__ABC01"),
        rank: Rank::Genus,
        action: Some(Action::Novel),
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
    }];
    assert_eq!(
        to_tsv(&rows),
        "taxon\trank\taction\tflag\tn_genomes_before\tn_genomes_after\tbootstrap_before\t\
         bootstrap_after\tred_before\tred_after\tprior_taxon\tnested_in\tsister_taxa\n\
         g__ABC01\tg\tnovel\tnovel_low_red,displaced_from_prior\t\t2\t\t100.0\t\t0.840\t\
         g__Old,g__Older\t\t\n"
    );
}
