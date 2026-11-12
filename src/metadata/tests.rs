use super::*;
use std::io::Cursor;

const TAX_A: &str = "d__Bacteria;p__Pseudomonadota;g__Escherichia";
const TAX_B: &str = "d__Bacteria;p__Bacillota;g__Bacillus";

fn parse(contents: &str) -> CuratomaticResult<Metadata> {
    Metadata::from_reader(Cursor::new(contents))
}

fn gid(gid: &str) -> GenomeId {
    GenomeId::new(gid)
}

fn ids(gids: &[&str]) -> BTreeSet<GenomeId> {
    gids.iter().map(|&gid| GenomeId::new(gid)).collect()
}

fn sample() -> Metadata {
    parse(&format!(
        "accession\tncbi_wgs_formatted\tgtdb_taxonomy\n\
         GB_GCA_000000001.1\tABC\t{TAX_A}\n\
         RS_GCF_000000002.1\tnone\t{TAX_B}\n\
         GCA_000000003.1\tDEF\tnone\n"
    ))
    .unwrap()
}

#[test]
fn parses_rows() {
    let meta = sample();
    assert_eq!(meta.rows.len(), 3);

    let row = &meta.rows[&gid("G000000001")];
    assert_eq!(row.accession, gid("G000000001"));
    assert_eq!(row.wgs.as_deref(), Some("ABC"));
    assert_eq!(
        row.taxonomy.as_ref().unwrap(),
        &Taxonomy::from_string(TAX_A).unwrap()
    );

    assert_eq!(meta.rows[&gid("G000000002")].wgs, None);
    assert!(meta.rows[&gid("G000000002")].taxonomy.is_some());
    assert_eq!(meta.rows[&gid("G000000003")].taxonomy, None);
}

#[test]
fn optional_columns_may_be_absent() {
    let meta = parse("accession\nGCA_000000001.1\n").unwrap();
    assert_eq!(meta.rows[&gid("G000000001")].taxonomy, None);
    assert_eq!(meta.rows[&gid("G000000001")].wgs, None);
}

#[test]
fn column_order_and_blank_lines() {
    let meta = parse(&format!(
        "gtdb_taxonomy\taccession\n{TAX_A}\tGCA_000000001.1\n\n"
    ))
    .unwrap();
    assert_eq!(meta.rows.len(), 1);
    assert!(meta.rows[&gid("G000000001")].taxonomy.is_some());
}

#[test]
fn missing_accession_column() {
    assert_eq!(
        parse("gtdb_taxonomy\n").unwrap_err(),
        CuratomaticError::MissingColumn("accession".to_owned())
    );
    assert!(parse("").is_err());
}

#[test]
fn short_row_is_an_error() {
    assert_eq!(
        parse("accession\tgtdb_taxonomy\nGCA_000000001.1\tnone\nGCA_000000002.1\tx\n")
            .map(|m| m.rows.len()),
        Err(CuratomaticError::UnknownRank("x".to_owned()))
    );
    assert_eq!(
        parse("accession\tncbi_wgs_formatted\tgtdb_taxonomy\nGCA_000000001.1\tA\n"),
        Err(CuratomaticError::MissingValue {
            line: 2,
            column: "gtdb_taxonomy".to_owned()
        })
    );
}

#[test]
fn missing_file() {
    let err = Metadata::from_file(Path::new("does/not/exist.tsv")).unwrap_err();
    assert!(matches!(err, CuratomaticError::Io(_)));
}

#[test]
fn subset() {
    let mut meta = sample();
    meta.subset_to_ids(&ids(&["G000000001", "G000000003", "G999"]));
    assert_eq!(
        meta.rows.keys().map(GenomeId::as_str).collect::<Vec<_>>(),
        ["G000000001", "G000000003"]
    );
}

#[test]
fn taxonomy_from_gids() {
    let meta = sample();
    let out = meta.get_taxonomy_from_gids(&ids(&["G000000001", "G000000003", "G999"]));
    // G000000003 has no taxonomy, and G999 is unknown.
    assert_eq!(out.keys().map(|gid| gid.as_str()).collect::<Vec<_>>(),
        ["G000000001"]);
}

#[test]
fn a_row_with_only_an_accession_skips_the_optional_columns() {
    // Python only reads the optional columns when `len(cols) > 1`.
    let meta = parse("accession\tncbi_wgs_formatted\tgtdb_taxonomy\nGCA_000000001.1\n").unwrap();
    assert_eq!(meta.rows[&gid("G000000001")].wgs, None);
    assert_eq!(meta.rows[&gid("G000000001")].taxonomy, None);
}

#[test]
fn none_is_only_a_placeholder_when_it_is_the_whole_value() {
    let meta = parse("accession\tncbi_wgs_formatted\nGCA_000000001.1\tNone\nGCA_000000002.1\tnone\n").unwrap();
    assert_eq!(meta.rows[&gid("G000000001")].wgs.as_deref(), Some("None"));
    assert_eq!(meta.rows[&gid("G000000002")].wgs, None);
}

#[test]
fn a_taxonomy_column_with_another_name_is_ignored() {
    // Python only reads `gtdb_taxonomy`.
    let meta = parse(&format!("accession\ttaxonomy\nGCA_000000001.1\t{TAX_A}\n")).unwrap();
    assert_eq!(meta.rows[&gid("G000000001")].taxonomy, None);
}

#[test]
fn the_last_row_wins_for_a_duplicate_accession() {
    let meta = parse("accession\tncbi_wgs_formatted\nGCA_000000001.1\tAAA\nGB_GCA_000000001.2\tBBB\n").unwrap();
    assert_eq!(meta.rows.len(), 1);
    assert_eq!(meta.rows[&gid("G000000001")].wgs.as_deref(), Some("BBB"));
}

#[test]
fn the_last_column_wins_for_a_duplicate_header_name() {
    let meta = parse("accession\tncbi_wgs_formatted\tncbi_wgs_formatted\nGCA_000000001.1\tAAA\tBBB\n").unwrap();
    assert_eq!(meta.rows[&gid("G000000001")].wgs.as_deref(), Some("BBB"));
}

#[test]
fn a_row_missing_the_wgs_value_names_the_column_and_line() {
    assert_eq!(
        parse("accession\tgtdb_taxonomy\tncbi_wgs_formatted\nGCA_000000001.1\tnone\n"),
        Err(CuratomaticError::MissingValue {
            line: 2,
            column: "ncbi_wgs_formatted".to_owned()
        })
    );
}

#[test]
fn a_row_too_short_for_the_accession_column_is_an_error() {
    assert_eq!(
        parse("gtdb_taxonomy\taccession\nnone\n"),
        Err(CuratomaticError::MissingValue {
            line: 2,
            column: "accession".to_owned()
        })
    );
}

#[test]
fn an_invalid_taxonomy_is_an_error() {
    assert_eq!(
        parse("accession\tgtdb_taxonomy\nGCA_000000001.1\tx__Foo\n"),
        Err(CuratomaticError::UnknownRank("x__Foo".to_owned()))
    );
}

#[test]
fn from_file_reads_a_file() {
    let path = std::env::temp_dir().join(format!("curatomatic_meta_{}.tsv", std::process::id()));
    std::fs::write(&path, format!("accession\tgtdb_taxonomy\nGCA_000000001.1\t{TAX_A}\n")).unwrap();
    let meta = Metadata::from_file(&path);
    std::fs::remove_file(&path).unwrap();
    assert!(meta.unwrap().rows[&gid("G000000001")].taxonomy.is_some());
}

#[test]
fn subset_to_no_ids_empties_the_metadata() {
    let mut meta = sample();
    meta.subset_to_ids(&BTreeSet::new());
    assert!(meta.rows.is_empty());
}

#[test]
fn an_empty_wgs_value_is_an_error() {
    for wgs in ["", " "] {
        assert_eq!(
            parse(&format!(
                "accession\tncbi_wgs_formatted\tgtdb_taxonomy\nGCA_000000001.1\tABC\tnone\nGCA_000000002.1\t{wgs}\tnone\n"
            )),
            Err(CuratomaticError::EmptyValue {
                line: 3,
                column: "ncbi_wgs_formatted".to_owned(),
            }),
            "wgs = {wgs:?}"
        );
    }
}
