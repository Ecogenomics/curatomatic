use super::*;

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

use std::path::Path;

use crate::genome::GenomeId;
use crate::taxonomy::Taxonomy;

fn metadata(contents: &str) -> Metadata {
    Metadata::from_reader(io::Cursor::new(contents), Path::new("names.tsv"), "wgs").unwrap()
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
            (Level::Info, "0 of 1,500 tree genomes have a prior taxonomy.".to_owned()),
        ]
    );
}

#[test]
fn coverage_warns_when_no_tree_genome_has_a_taxonomy() {
    let mut meta = metadata("accession\twgs\nG1\tnone\n");
    meta.set_taxonomy(GenomeId::new("G9"), Taxonomy::from_string("d__Bacteria").unwrap());
    meta.subset_to_ids(&[GenomeId::new("G1")].into());
    assert_eq!(
        taxonomy_coverage_messages(&meta, 2),
        vec![
            (Level::Info, "0 of 2 tree genomes have a prior taxonomy.".to_owned()),
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
    meta.set_taxonomy(GenomeId::new("G1"), Taxonomy::from_string("d__Bacteria").unwrap());
    assert_eq!(
        taxonomy_coverage_messages(&meta, 2),
        vec![(Level::Info, "1 of 2 tree genomes have a prior taxonomy.".to_owned())]
    );
}

fn parse_cli(args: &[&str]) -> Result<Cli, clap::Error> {
    Cli::try_parse_from(["curatomatic", "tree.nwk", "red.tsv", "out"].iter().chain(args))
}

#[test]
fn the_novel_name_file_and_column_are_given_together() {
    let cli = parse_cli(&[
        "--novel-name-file",
        "meta.tsv",
        "--novel-name-column",
        "ncbi_wgs_formatted",
    ])
    .unwrap();
    assert_eq!(cli.novel_name_file, Some(PathBuf::from("meta.tsv")));
    assert_eq!(cli.novel_name_column.as_deref(), Some("ncbi_wgs_formatted"));

    let cli = parse_cli(&[]).unwrap();
    assert_eq!((cli.novel_name_file, cli.novel_name_column), (None, None));
}

#[test]
fn the_novel_name_file_requires_a_column_and_the_column_requires_a_file() {
    for args in [
        &["--novel-name-file", "meta.tsv"][..],
        &["--novel-name-column", "ncbi_wgs_formatted"][..],
    ] {
        let err = parse_cli(args).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::MissingRequiredArgument, "{args:?}");
    }
}

#[test]
fn meta_is_no_longer_an_option() {
    assert!(parse_cli(&["--meta", "meta.tsv"]).is_err());
}
