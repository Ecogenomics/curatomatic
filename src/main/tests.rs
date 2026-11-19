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

    let cli = Cli::try_parse_from([
        "curatomatic".as_ref(),
        tree.as_os_str(),
        red.as_os_str(),
        out_dir.as_os_str(),
        "--taxonomy".as_ref(),
        path_taxonomy.as_os_str(),
    ])
    .unwrap();
    let result = curate(&cli);
    (out_dir, result)
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
fn taxa_check_accepts_an_unresolvable_novel_taxon_on_more_than_one_node() {
    let before = counts(&["d__D"]);
    let after = counts(&["d__D", "g__UNRESOLVABLE", "g__UNRESOLVABLE"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&["g__UNRESOLVABLE"])),
        Ok(())
    );
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
    let before = counts(&["d__D", "g__UNRESOLVABLE", "g__UNRESOLVABLE"]);
    let after = counts(&["d__D", "g__UNRESOLVABLE"]);
    assert_eq!(
        check_taxa_unchanged(&before, &after, &novel(&[])),
        taxa_changed("on a different number of nodes: [g__UNRESOLVABLE (2 -> 1)]")
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
