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

fn parse_cli(args: &[&str]) -> Result<Cli, clap::Error> {
    Cli::try_parse_from(
        ["curatomatic", "tree.nwk", "red.tsv", "out"]
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
fn the_taxonomy_file_is_optional() {
    let cli = parse_cli(&["--taxonomy", "taxonomy.tsv"]).unwrap();
    assert_eq!(cli.taxonomy, Some(PathBuf::from("taxonomy.tsv")));
    assert_eq!(parse_cli(&[]).unwrap().taxonomy, None);
}

/// Run `curate` on a small r237 style tree in a new directory, and return the output directory.
fn run_curate(name: &str, taxonomy: &str) -> (PathBuf, CuratomaticResult<()>) {
    let dir = std::env::temp_dir().join(format!("curatomatic-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let tree = dir.join("tree.nwk");
    std::fs::write(
        &tree,
        "((G1:0.1,G2:0.1)'100:g__A|red=0.95':0.1,G3:0.1)\
         'd__D;p__P;c__C;o__O;f__F|red=0.0000';\n",
    )
    .unwrap();
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
    // A leaf is labelled with its genome ID only, its taxa are in the taxonomy file.
    for gid in ["G1", "G2", "G3"] {
        assert!(tree.contains(&format!("{gid}:")), "{tree}");
    }
    assert!(!tree.contains("s__"), "{tree}");
    assert!(!tree.contains('|'), "{tree}");
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
