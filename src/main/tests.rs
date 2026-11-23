use super::*;
use crate::newick::NewickTree;

// Python formats counts with `f'{n:,}'`.
#[test]
fn separate_with_commas_matches_python_comma_formatting() {
    let cases = [
        (0_usize, "0"),
        (7, "7"),
        (999, "999"),
        (1000, "1,000"),
        (12_345, "12,345"),
        (100_000, "100,000"),
        (1_234_567, "1,234,567"),
    ];
    for (n, expected) in cases {
        assert_eq!(n.separate_with_commas(), expected, "n = {n}");
    }
}

use std::path::{Path, PathBuf};

use clap::Parser;

use crate::genome::GenomeId;
use crate::metadata::COL_ACCESSION;
use crate::taxonomy::Taxonomy;

fn metadata(contents: &str) -> Metadata {
    Metadata::from_reader(
        io::Cursor::new(contents),
        Path::new("names.tsv"),
        COL_ACCESSION,
        "wgs",
    )
    .unwrap()
}

#[test]
fn coverage_warns_when_no_prior_taxonomy_is_given() {
    let meta = metadata("accession\twgs\nG1\tABC\n");
    assert_eq!(
        taxonomy_coverage_messages(&meta, 1_500),
        vec![
            (
                Level::Warn,
                "No prior taxonomy was given, downward refinement is disabled.".to_owned()
            ),
            (
                Level::Info,
                "0 of 1,500 tree genomes have a prior taxonomy.".to_owned()
            ),
        ]
    );
}

#[test]
fn coverage_warns_when_no_tree_genome_has_a_taxonomy() {
    let mut meta = metadata("accession\twgs\nG1\tnone\n");
    meta.set_taxonomy(
        GenomeId::new("G9"),
        Taxonomy::from_string("d__Bacteria").unwrap(),
    );
    meta.subset_to_ids(&[GenomeId::new("G1")].into());
    assert_eq!(
        taxonomy_coverage_messages(&meta, 2),
        vec![
            (
                Level::Info,
                "0 of 2 tree genomes have a prior taxonomy.".to_owned()
            ),
            (
                Level::Warn,
                "No tree genomes have a prior taxonomy, downward refinement is disabled."
                    .to_owned()
            ),
        ]
    );
}

#[test]
fn coverage_only_reports_the_count_when_genomes_have_a_taxonomy() {
    let mut meta = metadata("accession\twgs\nG1\tABC\nG2\tnone\n");
    meta.set_taxonomy(
        GenomeId::new("G1"),
        Taxonomy::from_string("d__Bacteria").unwrap(),
    );
    assert_eq!(
        taxonomy_coverage_messages(&meta, 2),
        vec![(
            Level::Info,
            "1 of 2 tree genomes have a prior taxonomy.".to_owned()
        )]
    );
}

/// Run `curate` on a small r237 style tree in a new directory, and return the output directory.
fn run_curate(name: &str, taxonomy: &str) -> (PathBuf, CuratomaticResult<()>) {
    run_curate_on(
        name,
        "((G1:0.1,G2:0.1)'100:g__A|red=0.95':0.1,G3:0.1)\
         'd__D;p__P;c__C;o__O;f__F|red=0.0000';\n",
        taxonomy,
    )
}

/// Run `curate` on `newick` in a new directory, and return the output directory.
fn run_curate_on(name: &str, newick: &str, taxonomy: &str) -> (PathBuf, CuratomaticResult<()>) {
    run_curate_with(name, newick, taxonomy, &[])
}

/// Run `curate` on `newick` with `extra` arguments, see [`run_curate_on`].
fn run_curate_with(
    name: &str,
    newick: &str,
    taxonomy: &str,
    extra: &[&str],
) -> (PathBuf, CuratomaticResult<()>) {
    let dir = std::env::temp_dir().join(format!("curatomatic-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let tree = dir.join("tree.nwk");
    std::fs::write(&tree, newick).unwrap();
    let red = dir.join("red.tsv");
    std::fs::write(
        &red,
        "rank\tmedian\nPhylum\t0.3\nClass\t0.4\nOrder\t0.5\nFamily\t0.7\nGenus\t0.9\n",
    )
    .unwrap();
    let path_taxonomy = dir.join("taxonomy.tsv");
    std::fs::write(&path_taxonomy, taxonomy).unwrap();
    let out_dir = dir.join("out");
    std::fs::create_dir_all(&out_dir).unwrap();

    let mut args: Vec<&std::ffi::OsStr> = vec![
        "curatomatic".as_ref(),
        tree.as_os_str(),
        red.as_os_str(),
        out_dir.as_os_str(),
        "--taxonomy".as_ref(),
        path_taxonomy.as_os_str(),
    ];
    args.extend(extra.iter().map(std::ffi::OsStr::new));
    let cli = Cli::try_parse_from(args).unwrap();
    let result = curate(&cli);
    (out_dir, result)
}

/// `g__A` is on an unsupported node (bootstrap 50) above a supported clade of G1 and G2, and G3.
/// The decoration (the taxonomy file) puts every genome in `g__A`.
const MOVE_DOWN_TREE: &str = "(((G1:0.1,G2:0.1)'100|red=0.97':0.1,G3:0.1)'50:g__A|red=0.95':0.1)\
     'd__D;p__P;c__C;o__O;f__F|red=0.0000';\n";
const MOVE_DOWN_TAXONOMY: &str = "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
     G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two\n\
     G3\td__D;p__P;c__C;o__O;f__F;g__A;s__A three\n";

/// Run `curate` on [`MOVE_DOWN_TREE`], with `prior` as `--prior-taxonomy` if given, and return the
/// genus of each genome in the output.
fn genera_after_move_down(name: &str, prior: Option<&str>) -> Vec<String> {
    let prior_path = std::env::temp_dir().join(format!(
        "curatomatic-{name}-{}-prior.tsv",
        std::process::id()
    ));
    let prior_arg = prior_path.to_str().unwrap().to_owned();
    let mut extra = Vec::new();
    if let Some(prior) = prior {
        std::fs::write(&prior_path, prior).unwrap();
        extra = vec!["--prior-taxonomy", prior_arg.as_str()];
    }
    let (out_dir, result) = run_curate_with(name, MOVE_DOWN_TREE, MOVE_DOWN_TAXONOMY, &extra);
    result.unwrap();
    std::fs::read_to_string(out_dir.join("tree.curatomatic.taxonomy.tsv"))
        .unwrap()
        .lines()
        .map(|row| {
            let (gid, taxonomy) = row.split_once('\t').unwrap();
            let genus = taxonomy.split(';').find(|t| t.starts_with("g__")).unwrap();
            format!("{gid} {genus}")
        })
        .collect()
}

#[test]
fn a_taxon_is_not_moved_down_without_a_prior_taxonomy() {
    // The taxonomy file is the decoration of the same tree, so it can never move a taxon down.
    assert_eq!(
        genera_after_move_down("no-prior", None),
        ["G1 g__A", "G2 g__A", "G3 g__A"]
    );
}

#[test]
fn a_taxon_is_moved_down_to_the_genomes_of_the_prior_taxonomy() {
    let genera = genera_after_move_down(
        "prior",
        Some(
            "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
             G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two\n\
             G3\td__D;p__P;c__C;o__O;f__F;g__B;s__B three\n",
        ),
    );
    // G3 is left without a genus, and gets a novel one.
    assert_eq!(genera, ["G1 g__A", "G2 g__A", "G3 g__G3"]);
}

#[test]
fn a_genome_missing_from_the_prior_taxonomy_has_no_prior() {
    // G3 is not in the prior file: it is not an error, and G3 says nothing about g__A (its
    // species `s__A three` does not count, as `g__A` is not a Latin name, see rule D2a).
    let genera = genera_after_move_down(
        "prior-missing",
        Some(
            "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
             G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two\n\
             G9\td__D;p__P;c__C;o__O;f__F;g__A;s__A nine\n",
        ),
    );
    assert_eq!(genera, ["G1 g__A", "G2 g__A", "G3 g__G3"]);
}

#[test]
fn curate_writes_one_seven_rank_row_per_genome_from_the_final_tree() {
    let (out_dir, result) = run_curate(
        "writes",
        "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;s__{unresolved} sp2\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__B;s__B three\n\
         G9\td__D;p__P;c__C;o__O;f__F;g__B;s__B nine\n",
    );
    result.unwrap();
    let taxonomy = std::fs::read_to_string(out_dir.join("tree.curatomatic.taxonomy.tsv")).unwrap();
    assert_eq!(
        taxonomy,
        "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;s__{unresolved} sp2\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__B;s__B three\n"
    );
    let tree = std::fs::read_to_string(out_dir.join("tree.curatomatic.tree")).unwrap();
    // Leaf labels are never changed: the taxa of a leaf are in the taxonomy file only.
    let leaf_labels = |text: &str| {
        let newick = NewickTree::parse(text).unwrap();
        newick
            .leaves()
            .map(|leaf| newick.label(leaf).unwrap().to_owned())
            .collect::<Vec<_>>()
    };
    let tree_in = std::fs::read_to_string(out_dir.join("../tree.nwk")).unwrap();
    assert_eq!(leaf_labels(&tree), leaf_labels(&tree_in));
    assert_eq!(leaf_labels(&tree), ["G1", "G2", "G3"]);
}

#[test]
fn curate_fails_when_a_tree_genome_is_not_in_the_taxonomy_file() {
    let (_, result) = run_curate("missing", "G1\td__D;s__A one\n");
    assert_eq!(
        result,
        Err(CuratomaticError::MissingFromTaxonomyFile {
            n: 2,
            first: GenomeId::new("G2")
        })
    );
}

#[test]
fn curate_keeps_a_single_genome_genus_nested_in_a_tree_genus() {
    // G2's genus is only in the file (it was on its dummy node), inside the tree's `g__A`.
    let (out_dir, result) = run_curate(
        "nested",
        "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;g__N;s__N two\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__B;s__B three\n",
    );
    result.unwrap();
    let taxonomy = std::fs::read_to_string(out_dir.join("tree.curatomatic.taxonomy.tsv")).unwrap();
    assert!(
        taxonomy.contains("G2\td__D;p__P;c__C;o__O;f__F;g__N;s__N two\n"),
        "{taxonomy}"
    );
}

#[test]
fn curate_writes_the_taxa_held_by_each_leaf_to_the_leaf_taxonomy_file() {
    // Only the taxa that the Newick tree can't hold: no placeholders and no taxa on the path.
    let (out_dir, result) = run_curate(
        "leaf-taxonomy",
        "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;g__N;s__N two\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__B;s__B three\n",
    );
    result.unwrap();
    let leaf_taxonomy =
        std::fs::read_to_string(out_dir.join("tree.curatomatic.leaf-taxonomy.tsv")).unwrap();
    assert_eq!(
        leaf_taxonomy,
        "G1\ts__A one\n\
         G2\tg__N;s__N two\n\
         G3\tg__B;s__B three\n"
    );
}

#[test]
fn curate_fails_when_the_taxonomy_file_disagrees_with_the_tree() {
    let (_, result) = run_curate(
        "disagrees",
        "G1\td__D;p__P;c__C;o__O;f__F;g__Other;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__B;s__B three\n",
    );
    assert_eq!(
        result,
        Err(CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G1"),
            message: "the tree places it in g__A, which is not in its row of the taxonomy file"
                .to_owned()
        })
    );
}

#[test]
fn curate_fails_without_writing_anything_when_the_file_disagrees_above_the_anchor() {
    // G1's nearest taxon on the tree (g__A) is in its row, but the tree puts it in p__P.
    let (out_dir, result) = run_curate(
        "disagrees-above-anchor",
        "G1\td__D;p__Q;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__B;s__B three\n",
    );
    assert_eq!(
        result,
        Err(CuratomaticError::TaxonomyFileConflict {
            gid: GenomeId::new("G1"),
            message: "the tree places it in p__P, which is not in its row of the taxonomy file"
                .to_owned()
        })
    );
    assert_eq!(std::fs::read_dir(out_dir).unwrap().count(), 0);
}

#[test]
fn curate_fails_without_writing_anything_when_a_leaf_has_no_domain() {
    let (out_dir, result) = run_curate_on(
        "no-domain",
        "((G1:0.1,G2:0.1)'100:g__A|red=0.95':0.1,G3:0.1)'p__P;c__C;o__O;f__F|red=0.0000';\n",
        "G1\td__;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__;p__P;c__C;o__O;f__F;g__A;s__A two\n\
         G3\td__;p__P;c__C;o__O;f__F;g__B;s__B three\n",
    );
    assert_eq!(
        result,
        Err(CuratomaticError::LeafWithoutDomain {
            n: 3,
            first: GenomeId::new("G1"),
        })
    );
    assert_eq!(std::fs::read_dir(out_dir).unwrap().count(), 0);
}

#[test]
fn curate_fails_without_writing_anything_when_a_leaf_has_no_species() {
    let (out_dir, result) = run_curate(
        "no-species",
        "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;s__\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__B;s__\n",
    );
    assert_eq!(
        result,
        Err(CuratomaticError::LeafWithoutSpecies {
            n: 2,
            first: GenomeId::new("G2"),
        })
    );
    assert_eq!(std::fs::read_dir(out_dir).unwrap().count(), 0);
}

const DEBUG_TREE: &str = "((G1:0.1,G2:0.1)'100:g__A|red=0.95':0.1,G3:0.1)\
                          'd__D;p__P;c__C;o__O;f__F|red=0.0000';\n";

/// A taxonomy file for [`DEBUG_TREE`], where G2 is in a genus nested in `g__A`, and G3 in a
/// family nested in `f__F`.
const DEBUG_TAXONOMY: &str = "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
                              G2\td__D;p__P;c__C;o__O;f__F;g__A;g__N;s__N two\n\
                              G3\td__D;p__P;c__C;o__O;f__F;f__Q;g__Q;s__Q three\n";

#[test]
fn curate_with_debug_writes_the_inputs_as_loaded_with_dummy_nodes() {
    let (out_dir, result) = run_curate_with("debug", DEBUG_TREE, DEBUG_TAXONOMY, &["--debug"]);
    result.unwrap();
    let zombie = std::fs::read_to_string(out_dir.join("tree.debug.zombie.tree")).unwrap();
    assert_eq!(
        zombie,
        "(((G1|RED=1.0:0.0,D-G1|RED=1.0:0.0)'s__A one|RED=1.0':0.1,\
         (G2|RED=1.0:0.0,D-G2|RED=1.0:0.0)'g__N; s__N two|RED=1.0':0.1)'100:g__A|red=0.95':0.1,\
         (G3|RED=1.0:0.0,D-G3|RED=1.0:0.0)'f__Q; g__Q; s__Q three|RED=1.0':0.1)\
         'd__D;p__P;c__C;o__O;f__F|red=0.0000';\n"
    );
    // The taxonomy as loaded is the input file, before any curation.
    let taxonomy = std::fs::read_to_string(out_dir.join("tree.debug.taxonomy.tsv")).unwrap();
    assert_eq!(taxonomy, DEBUG_TAXONOMY);
}

#[test]
fn curate_without_debug_writes_no_debug_files() {
    let (out_dir, result) = run_curate_on("no-debug", DEBUG_TREE, DEBUG_TAXONOMY);
    result.unwrap();
    for name in ["tree.debug.zombie.tree", "tree.debug.taxonomy.tsv"] {
        assert!(!out_dir.join(name).exists(), "{name}");
    }
}

// The rule: curatomatic only ever moves taxa or creates novel taxa. It never removes a taxon, and
// never edits the text of one (see `AGENTS.md`).

fn counts(taxa: &[&str]) -> TaxonCounts {
    let taxa: Vec<Taxon> = taxa.iter().map(|&taxon| Taxon::new(taxon)).collect();
    TaxonCounts::from_taxa(&taxa)
}

fn novel(taxa: &[&str]) -> BTreeSet<Taxon> {
    taxa.iter().map(|&taxon| Taxon::new(taxon)).collect()
}

fn taxa_changed(message: &str) -> CuratomaticResult<()> {
    Err(CuratomaticError::TaxaChanged(message.to_owned()))
}

#[test]
fn taxa_check_accepts_taxa_that_were_only_moved() {
    // A move changes the node of a taxon, never how many nodes hold it.
    let taxa = counts(&["d__D", "f__F", "g__A", "s__A one"]);
    assert_eq!(check_taxa_unchanged(&taxa, &taxa, &novel(&[])), Ok(()));
}

#[test]
fn taxa_check_accepts_novel_taxa() {
    let before = counts(&["d__D", "g__A"]);
    let after = counts(&["d__D", "g__A", "f__G1", "g__G2"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&["f__G1", "g__G2"])),
        Ok(())
    );
}

#[test]
fn taxa_check_fails_when_a_novel_taxon_is_on_more_than_one_node() {
    // The unresolvable name is never placed, so it has no exception.
    for name in ["g__G1", "g__UNRESOLVABLE"] {
        let before = counts(&["d__D"]);
        let after = counts(&["d__D", name, name]);
        assert_eq!(
            check_taxa_unchanged(&before, &after, &novel(&[name])),
            taxa_changed(&format!("novel, but on more than one node: [{name} (2)]"))
        );
    }
}

#[test]
fn taxa_check_ignores_placeholders() {
    // A placeholder is not a name, and can be refined freely.
    let before = counts(&["d__D", "p__", "g__", "g__A"]);
    let after = counts(&["d__D", "g__A", "p__G1"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&["p__G1"])),
        Ok(())
    );
}

#[test]
fn taxa_check_fails_when_a_taxon_is_removed() {
    let before = counts(&["d__D", "f__F", "g__A"]);
    let after = counts(&["d__D", "g__A"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&[])),
        taxa_changed("removed: [f__F]")
    );
}

#[test]
fn taxa_check_fails_when_the_text_of_a_taxon_is_edited() {
    // An edited taxon is the old text removed, and new text that is not a novel taxon.
    let before = counts(&["d__D", "g__Escherichia"]);
    let after = counts(&["d__D", "g__Escherichia_A"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&[])),
        taxa_changed(
            "removed: [g__Escherichia]; gained, but not created as a novel taxon: \
             [g__Escherichia_A]"
        )
    );
}

#[test]
fn taxa_check_fails_when_an_unresolved_species_is_edited() {
    let before = counts(&["d__D", "s__{unresolved} sp964401205"]);
    let after = counts(&["d__D", "s__{unresolved} G964401205"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&[])),
        taxa_changed(
            "removed: [s__{unresolved} sp964401205]; gained, but not created as a novel taxon: \
             [s__{unresolved} G964401205]"
        )
    );
}

#[test]
fn taxa_check_fails_when_a_moved_taxon_is_left_behind() {
    // Copied instead of moved: it is on both the old and the new node.
    let before = counts(&["d__D", "g__A"]);
    let after = counts(&["d__D", "g__A", "g__A"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&[])),
        taxa_changed("on a different number of nodes: [g__A (1 -> 2)]")
    );
}

#[test]
fn taxa_check_fails_when_one_copy_of_a_taxon_is_removed() {
    let before = counts(&["d__D", "g__A", "g__A"]);
    let after = counts(&["d__D", "g__A"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&[])),
        taxa_changed("on a different number of nodes: [g__A (2 -> 1)]")
    );
}

#[test]
fn taxa_check_fails_when_a_taxon_is_added_that_is_not_novel() {
    let before = counts(&["d__D", "g__A"]);
    let after = counts(&["d__D", "g__A", "g__B"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&[])),
        taxa_changed("gained, but not created as a novel taxon: [g__B]")
    );
}

#[test]
fn taxa_check_fails_when_a_novel_taxon_reuses_an_existing_name() {
    let before = counts(&["d__D", "g__A"]);
    let after = counts(&["d__D", "g__A", "g__A"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&["g__A"])),
        taxa_changed(
            "on a different number of nodes: [g__A (1 -> 2)]; reused by a novel taxon: [g__A]"
        )
    );
}

#[test]
fn taxa_check_lists_every_problem() {
    let before = counts(&["d__D", "f__F", "g__A", "g__B"]);
    let after = counts(&["d__D", "g__A", "g__A", "g__C", "g__N"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&["g__N", "g__B"])),
        taxa_changed(
            "removed: [f__F, g__B]; on a different number of nodes: [g__A (1 -> 2)]; \
             reused by a novel taxon: [g__B]; gained, but not created as a novel taxon: [g__C]"
        )
    );
}

#[test]
fn taxa_changed_error_states_the_rule() {
    let error = CuratomaticError::TaxaChanged("removed: [f__F]".to_owned());
    assert_eq!(
        error.to_string(),
        "Curatomatic may only move taxa or create novel taxa, never remove them or edit their \
         text, but the taxa on the tree changed: removed: [f__F]"
    );
}

#[test]
fn curate_passes_the_taxa_check_after_moving_a_taxon_and_creating_one() {
    // g__A is moved up (RED 0.97 to 0.92, closer to the genus median of 0.9), taking in G3, and a
    // novel genus is created for G4. The run, and so the check at its end, must succeed.
    let (out_dir, result) = run_curate_on(
        "taxa-check",
        "(((G1:0.1,G2:0.1)'100:g__A|red=0.97':0.1,G3:0.1)'100|red=0.92':0.1,G4:0.1)\
         'd__D;p__P;c__C;o__O;f__F|red=0.0000';\n",
        "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__;s__{unresolved} sp3\n\
         G4\td__D;p__P;c__C;o__O;f__F;g__;s__{unresolved} sp4\n",
    );
    result.unwrap();
    let taxonomy = std::fs::read_to_string(out_dir.join("tree.curatomatic.taxonomy.tsv")).unwrap();
    assert_eq!(
        taxonomy,
        "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
         G2\td__D;p__P;c__C;o__O;f__F;g__A;s__A two\n\
         G3\td__D;p__P;c__C;o__O;f__F;g__A;s__{unresolved} sp3\n\
         G4\td__D;p__P;c__C;o__O;f__F;g__G4;s__{unresolved} sp4\n"
    );
}

// Any input may be gzipped (its name ends in `.gz`): the run must be the same as on the plain
// files, outputs and their names included.

/// The inputs of the README example, in `example/`: `--taxonomy`, `--prior-taxonomy`,
/// `--novel-name-file`, `TREE` and `MEDIAN_RED`.
const EXAMPLE_INPUTS: [&str; 5] = [
    "gtdb_r237_decorated_taxonomy.tsv",
    "gtdb_r232_bac120_taxonomy.tsv",
    "gtdb_r237_bac120_metadata_slim.tsv",
    "gtdb_r237_median_red_marked.nwk",
    "gtdb_r237_median_red_rank_distribution.tsv",
];

/// The arguments of the README example, on the inputs `paths` (in the order of
/// [`EXAMPLE_INPUTS`]), writing to `out_dir`.
fn example_cli(paths: &[PathBuf], out_dir: &Path) -> Cli {
    let [taxonomy, prior, names, tree, red] = paths else {
        panic!("expected the five example inputs");
    };
    Cli::try_parse_from([
        "curatomatic".as_ref(),
        "--min-bs".as_ref(),
        "95.0".as_ref(),
        "--taxonomy".as_ref(),
        taxonomy.as_os_str(),
        "--prior-taxonomy".as_ref(),
        prior.as_os_str(),
        "--novel-name-file".as_ref(),
        names.as_os_str(),
        "--novel-accession".as_ref(),
        "db_name".as_ref(),
        "--novel-name-column".as_ref(),
        "ncbi_wgs_formatted".as_ref(),
        tree.as_os_str(),
        red.as_os_str(),
        out_dir.as_os_str(),
    ])
    .unwrap()
}

/// The name and contents of every file in `dir`, ordered by name.
fn files_in(dir: &Path) -> std::collections::BTreeMap<String, Vec<u8>> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_string_lossy().into_owned();
            (name, std::fs::read(&path).unwrap())
        })
        .collect()
}

#[test]
fn curate_on_the_gzipped_example_inputs_writes_the_outputs_of_the_plain_inputs() {
    use flate2::Compression;
    use flate2::write::GzEncoder;

    let example = Path::new(env!("CARGO_MANIFEST_DIR")).join("example");
    let dir = std::env::temp_dir().join(format!(
        "curatomatic-gzipped-example-{}",
        std::process::id()
    ));
    if dir.exists() {
        std::fs::remove_dir_all(&dir).unwrap();
    }
    std::fs::create_dir_all(&dir).unwrap();

    let plain: Vec<PathBuf> = EXAMPLE_INPUTS
        .iter()
        .map(|name| example.join(name))
        .collect();
    let gzipped: Vec<PathBuf> = EXAMPLE_INPUTS
        .iter()
        .map(|name| {
            let path = dir.join(format!("{name}.gz"));
            let mut encoder =
                GzEncoder::new(std::fs::File::create(&path).unwrap(), Compression::fast());
            io::copy(
                &mut std::fs::File::open(example.join(name)).unwrap(),
                &mut encoder,
            )
            .unwrap();
            encoder.finish().unwrap();
            path
        })
        .collect();

    // Both runs at once: each takes a while on the full example.
    let (out_plain, out_gzipped) = (dir.join("out-plain"), dir.join("out-gzipped"));
    std::fs::create_dir_all(&out_plain).unwrap();
    std::fs::create_dir_all(&out_gzipped).unwrap();
    let (result_plain, result_gzipped) = std::thread::scope(|scope| {
        let plain_run = scope.spawn(|| curate(&example_cli(&plain, &out_plain)));
        let gzipped_run = scope.spawn(|| curate(&example_cli(&gzipped, &out_gzipped)));
        (plain_run.join().unwrap(), gzipped_run.join().unwrap())
    });
    result_plain.unwrap();
    result_gzipped.unwrap();

    let (files_plain, files_gzipped) = (files_in(&out_plain), files_in(&out_gzipped));
    std::fs::remove_dir_all(&dir).unwrap();
    // The outputs are named after the tree without `.gz`, and are not compressed.
    assert_eq!(
        files_plain.keys().collect::<Vec<_>>(),
        [
            "gtdb_r237_median_red_marked.curatomatic.leaf-taxonomy.tsv",
            "gtdb_r237_median_red_marked.curatomatic.report.tsv",
            "gtdb_r237_median_red_marked.curatomatic.taxonomy.tsv",
            "gtdb_r237_median_red_marked.curatomatic.tree",
        ]
    );
    assert_eq!(
        files_gzipped.keys().collect::<Vec<_>>(),
        files_plain.keys().collect::<Vec<_>>()
    );
    for (name, contents) in &files_plain {
        assert!(contents == &files_gzipped[name], "{name} differs");
    }
}
