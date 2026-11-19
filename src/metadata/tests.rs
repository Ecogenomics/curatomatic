use super::*;
use std::io::Cursor;

const TAX_A: &str = "d__Bacteria;p__Pseudomonadota;g__Escherichia";
const TAX_B: &str = "d__Bacteria;p__Bacillota;g__Bacillus";
const COL: &str = "ncbi_wgs_formatted";

fn parse_column(contents: &str, column: &str) -> CuratomaticResult<Metadata> {
    Metadata::from_reader(Cursor::new(contents), Path::new("names.tsv"), column)
}

fn parse(contents: &str) -> CuratomaticResult<Metadata> {
    parse_column(contents, COL)
}

fn gid(gid: &str) -> GenomeId {
    GenomeId::new(gid)
}

fn ids(gids: &[&str]) -> BTreeSet<GenomeId> {
    gids.iter().map(|&gid| GenomeId::new(gid)).collect()
}

fn wgs<'a>(meta: &'a Metadata, id: &str) -> Option<&'a str> {
    meta.rows[&gid(id)].wgs.as_deref()
}

fn sample() -> Metadata {
    let mut meta = parse(
        "accession\tncbi_wgs_formatted\tgtdb_taxonomy\n\
         GB_GCA_000000001.1\tABC\tignored\n\
         RS_GCF_000000002.1\tnone\tignored\n\
         GCA_000000003.1\tDEF\tignored\n",
    )
    .unwrap();
    meta.set_taxonomy(gid("G000000001"), Taxonomy::from_string(TAX_A).unwrap());
    meta.set_taxonomy(gid("G000000002"), Taxonomy::from_string(TAX_B).unwrap());
    meta
}

#[test]
fn parses_rows() {
    let meta =
        parse("accession\tncbi_wgs_formatted\nGB_GCA_000000001.1\tABC\nRS_GCF_000000002.1\tnone\n")
            .unwrap();
    assert_eq!(meta.rows.len(), 2);
    let row = &meta.rows[&gid("G000000001")];
    assert_eq!(row.accession, gid("G000000001"));
    assert_eq!(row.wgs.as_deref(), Some("ABC"));
    assert_eq!(row.taxonomy, None);
    assert_eq!(wgs(&meta, "G000000002"), None);
}

#[test]
fn the_file_never_gives_a_prior_taxonomy() {
    let meta = parse(&format!(
        "accession\tgtdb_taxonomy\tncbi_wgs_formatted\nGCA_000000001.1\t{TAX_A}\tABC\n"
    ))
    .unwrap();
    assert!(!meta.has_prior_taxonomy);
    assert_eq!(meta.n_with_taxonomy(), 0);
    assert_eq!(meta.rows[&gid("G000000001")].taxonomy, None);
}

#[test]
fn reads_a_custom_column_and_ignores_the_others() {
    let meta = parse_column(
        "extra\tmy_names\taccession\tncbi_wgs_formatted\n\
         x\tFOO\tGCA_000000001.1\tABC\n",
        "my_names",
    )
    .unwrap();
    assert_eq!(wgs(&meta, "G000000001"), Some("FOO"));
}

#[test]
fn a_missing_name_column_is_an_error_with_the_path() {
    let err = parse_column(
        "accession\tncbi_wgs_formatted\nGCA_000000001.1\tABC\n",
        "my_names",
    )
    .unwrap_err();
    assert_eq!(
        err,
        CuratomaticError::MissingColumn {
            path: "names.tsv".to_owned(),
            column: "my_names".to_owned()
        }
    );
    assert_eq!(
        err.to_string(),
        "Missing required column my_names in: names.tsv"
    );
}

#[test]
fn a_missing_accession_column_is_an_error_with_the_path() {
    let err = CuratomaticError::MissingColumn {
        path: "names.tsv".to_owned(),
        column: "accession".to_owned(),
    };
    assert_eq!(parse("ncbi_wgs_formatted\nABC\n").unwrap_err(), err);
    assert_eq!(parse("").unwrap_err(), err);
}

#[test]
fn each_placeholder_means_no_name() {
    for placeholder in [
        "none", "None", "NONE", "n/a", "N/A", "null", "NULL", "", " ", " none ",
    ] {
        let meta = parse(&format!(
            "accession\tncbi_wgs_formatted\tother\nGCA_000000001.1\t{placeholder}\tx\n"
        ))
        .unwrap();
        assert_eq!(meta.rows.len(), 1, "{placeholder:?}");
        assert_eq!(wgs(&meta, "G000000001"), None, "{placeholder:?}");
    }
}

#[test]
fn an_empty_last_value_is_no_name() {
    // The line is not trimmed as a whole, so the trailing empty value keeps its column.
    let meta =
        parse("accession\tncbi_wgs_formatted\nGCA_000000001.1\t\nGCA_000000002.1\t\r\n").unwrap();
    assert_eq!(wgs(&meta, "G000000001"), None);
    assert_eq!(wgs(&meta, "G000000002"), None);
}

#[test]
fn names_are_trimmed() {
    let meta = parse("accession\tncbi_wgs_formatted\n GCA_000000001.1 \t ABC \n").unwrap();
    assert_eq!(wgs(&meta, "G000000001"), Some("ABC"));
}

#[test]
fn a_placeholder_must_be_the_whole_value() {
    let meta = parse("accession\tncbi_wgs_formatted\nGCA_000000001.1\tnoneX\n").unwrap();
    assert_eq!(wgs(&meta, "G000000001"), Some("noneX"));
}

#[test]
fn column_order_and_blank_lines() {
    let meta = parse("ncbi_wgs_formatted\taccession\nABC\tGCA_000000001.1\n\n").unwrap();
    assert_eq!(meta.rows.len(), 1);
    assert_eq!(wgs(&meta, "G000000001"), Some("ABC"));
}

#[test]
fn a_row_missing_the_name_value_names_the_column_and_line() {
    assert_eq!(
        parse("accession\tother\tncbi_wgs_formatted\nGCA_000000001.1\tx\tA\nGCA_000000002.1\tx\n"),
        Err(CuratomaticError::MissingValue {
            line: 3,
            column: COL.to_owned()
        })
    );
}

#[test]
fn a_row_too_short_for_the_accession_column_is_an_error() {
    assert_eq!(
        parse("ncbi_wgs_formatted\taccession\nABC\n"),
        Err(CuratomaticError::MissingValue {
            line: 2,
            column: "accession".to_owned()
        })
    );
}

#[test]
fn missing_file() {
    let err = Metadata::from_file(Path::new("does/not/exist.tsv"), COL).unwrap_err();
    assert!(matches!(err, CuratomaticError::Io(message) if message.contains("does/not/exist.tsv")));
}

#[test]
fn from_file_reads_a_file() {
    let path = std::env::temp_dir().join(format!("curatomatic_meta_{}.tsv", std::process::id()));
    std::fs::write(&path, "accession\tmy_names\nGCA_000000001.1\tABC\n").unwrap();
    let meta = Metadata::from_file(&path, "my_names");
    std::fs::remove_file(&path).unwrap();
    assert_eq!(wgs(&meta.unwrap(), "G000000001"), Some("ABC"));
}

#[test]
fn the_last_row_wins_for_a_duplicate_accession() {
    let meta =
        parse("accession\tncbi_wgs_formatted\nGCA_000000001.1\tAAA\nGB_GCA_000000001.2\tBBB\n")
            .unwrap();
    assert_eq!(meta.rows.len(), 1);
    assert_eq!(wgs(&meta, "G000000001"), Some("BBB"));
}

#[test]
fn the_last_column_wins_for_a_duplicate_header_name() {
    let meta =
        parse("accession\tncbi_wgs_formatted\tncbi_wgs_formatted\nGCA_000000001.1\tAAA\tBBB\n")
            .unwrap();
    assert_eq!(wgs(&meta, "G000000001"), Some("BBB"));
}

#[test]
fn set_taxonomy_adds_a_row_and_marks_a_prior_taxonomy() {
    let mut meta = Metadata::default();
    assert!(!meta.has_prior_taxonomy);
    meta.set_taxonomy(gid("G1"), Taxonomy::from_string(TAX_A).unwrap());
    assert!(meta.has_prior_taxonomy);
    assert!(meta.has_taxonomy(&gid("G1")));
    assert_eq!(meta.rows[&gid("G1")].wgs, None);
}

#[test]
fn set_taxonomy_keeps_the_name() {
    let meta = sample();
    assert_eq!(wgs(&meta, "G000000001"), Some("ABC"));
    assert_eq!(
        meta.rows[&gid("G000000001")].taxonomy.as_ref().unwrap(),
        &Taxonomy::from_string(TAX_A).unwrap()
    );
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
fn subset_to_no_ids_empties_the_metadata() {
    let mut meta = sample();
    meta.subset_to_ids(&BTreeSet::new());
    assert!(meta.rows.is_empty());
}

#[test]
fn taxonomy_from_gids() {
    let meta = sample();
    let out = meta.get_taxonomy_from_gids(&ids(&["G000000001", "G000000003", "G999"]));
    // G000000003 has no taxonomy, and G999 is unknown.
    assert_eq!(
        out.keys().map(|gid| gid.as_str()).collect::<Vec<_>>(),
        ["G000000001"]
    );
}

#[test]
fn n_with_taxonomy_counts_rows_with_a_taxonomy() {
    let meta = sample();
    assert!(meta.has_prior_taxonomy);
    assert_eq!(meta.n_with_taxonomy(), 2);
}
