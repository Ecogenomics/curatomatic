use std::path::PathBuf;

use super::*;

fn parse_cli(args: &[&str]) -> Result<Cli, clap::Error> {
    Cli::try_parse_from(
        [
            "curatomatic",
            "tree.nwk",
            "red.tsv",
            "out",
            "--taxonomy",
            "taxonomy.tsv",
        ]
        .iter()
        .chain(args),
    )
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
        assert_eq!(
            err.kind(),
            clap::error::ErrorKind::MissingRequiredArgument,
            "{args:?}"
        );
    }
}

#[test]
fn the_novel_accession_column_defaults_to_accession() {
    assert_eq!(parse_cli(&[]).unwrap().novel_accession, "accession");
    let cli = parse_cli(&[
        "--novel-name-file",
        "meta.tsv",
        "--novel-name-column",
        "ncbi_wgs_formatted",
    ])
    .unwrap();
    assert_eq!(cli.novel_accession, "accession");
}

#[test]
fn the_novel_accession_column_is_parsed() {
    let cli = parse_cli(&[
        "--novel-name-file",
        "meta.tsv",
        "--novel-accession",
        "db_name",
        "--novel-name-column",
        "ncbi_wgs_formatted",
    ])
    .unwrap();
    assert_eq!(cli.novel_accession, "db_name");
}

#[test]
fn the_novel_accession_column_requires_a_file() {
    let err = parse_cli(&["--novel-accession", "db_name"]).unwrap_err();
    assert_eq!(err.kind(), clap::error::ErrorKind::MissingRequiredArgument);
}

#[test]
fn meta_is_no_longer_an_option() {
    assert!(parse_cli(&["--meta", "meta.tsv"]).is_err());
}

#[test]
fn the_taxonomy_file_is_required() {
    assert_eq!(
        parse_cli(&[]).unwrap().taxonomy,
        PathBuf::from("taxonomy.tsv")
    );
    // Without it, no leaf would hold its species.
    assert!(Cli::try_parse_from(["curatomatic", "tree.nwk", "red.tsv", "out"]).is_err());
}

#[test]
fn debug_is_off_unless_given() {
    assert!(!parse_cli(&[]).unwrap().debug);
    assert!(parse_cli(&["--debug"]).unwrap().debug);
}
