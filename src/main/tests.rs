use super::*;

// Python formats counts with `f'{n:,}'`.
#[test]
fn thousands_matches_python_comma_formatting() {
    let cases = [
        (0, "0"),
        (7, "7"),
        (999, "999"),
        (1000, "1,000"),
        (12_345, "12,345"),
        (100_000, "100,000"),
        (1_234_567, "1,234,567"),
    ];
    for (n, expected) in cases {
        assert_eq!(thousands(n), expected, "n = {n}");
    }
}

fn metadata(contents: &str) -> Metadata {
    Metadata::from_reader(io::Cursor::new(contents)).unwrap()
}

#[test]
fn coverage_warns_when_the_gtdb_taxonomy_column_is_missing() {
    let meta = metadata("accession\ttaxonomy\nG1\td__Bacteria\n");
    assert_eq!(
        taxonomy_coverage_messages(&meta, 1_500),
        vec![
            (
                Level::Warn,
                "The metadata file has no gtdb_taxonomy column, downward refinement is disabled."
                    .to_owned()
            ),
            (
                Level::Info,
                "0 of 1,500 tree genomes have a gtdb_taxonomy.".to_owned()
            ),
        ]
    );
}

#[test]
fn coverage_warns_when_no_genome_has_a_taxonomy() {
    let meta = metadata("accession\tgtdb_taxonomy\nG1\tnone\n");
    assert_eq!(
        taxonomy_coverage_messages(&meta, 2),
        vec![
            (Level::Info, "0 of 2 tree genomes have a gtdb_taxonomy.".to_owned()),
            (
                Level::Warn,
                "No tree genomes have a gtdb_taxonomy, downward refinement is disabled."
                    .to_owned()
            ),
        ]
    );
}

#[test]
fn coverage_only_reports_the_count_when_genomes_have_a_taxonomy() {
    let meta = metadata("accession\tgtdb_taxonomy\nG1\td__Bacteria\nG2\tnone\n");
    assert_eq!(
        taxonomy_coverage_messages(&meta, 2),
        vec![(Level::Info, "1 of 2 tree genomes have a gtdb_taxonomy.".to_owned())]
    );
}
