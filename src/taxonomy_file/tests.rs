use super::*;

fn read(contents: &str) -> CuratomaticResult<TaxonomyFile> {
    TaxonomyFile::from_reader(contents.as_bytes(), Path::new("taxonomy.tsv"))
}

const ECOLI: &str = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;\
                     f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli";

#[test]
fn reads_the_genome_and_its_taxonomy() {
    let file = read(&format!("G003697165\t{ECOLI}\n")).unwrap();
    let taxonomy = file.get(&GenomeId::new("G003697165")).unwrap();
    assert_eq!(taxonomy, &Taxonomy::from_string(ECOLI).unwrap());
}

#[test]
fn genome_ids_are_made_canonical_like_the_metadata() {
    let file = read(&format!(
        "GB_GCA_003697165.1\t{ECOLI}\nRS_GCF_000742135.1\td__Bacteria\n"
    ))
    .unwrap();
    let gids: Vec<&str> = file.iter().map(|(gid, _)| gid.as_str()).collect();
    assert_eq!(gids, ["G000742135", "G003697165"]);
}

#[test]
fn placeholders_and_unresolved_species_are_kept_as_read() {
    let file = read("G1\td__Bacteria;p__;g__;s__{unresolved} sp053939565\n").unwrap();
    let taxonomy = file.get(&GenomeId::new("G1")).unwrap();
    assert_eq!(taxonomy.p, vec![Taxon::new("p__")]);
    assert_eq!(taxonomy.s, vec![Taxon::new("s__{unresolved} sp053939565")]);
}

#[test]
fn nested_taxa_are_kept() {
    let file = read("G1\td__Bacteria;p__Bacillota;p__Bacillota_I;c__Bacilli\n").unwrap();
    assert_eq!(
        file.get(&GenomeId::new("G1")).unwrap().p,
        vec![Taxon::new("p__Bacillota"), Taxon::new("p__Bacillota_I")]
    );
}

#[test]
fn blank_lines_and_windows_line_endings_are_ignored() {
    let file = read("\nG1\td__Bacteria;s__A b\r\n\n").unwrap();
    assert_eq!(
        file.get(&GenomeId::new("G1")).unwrap().s,
        vec![Taxon::new("s__A b")]
    );
}

#[test]
fn a_line_without_two_columns_is_an_error() {
    for (contents, line) in [("G1\n", 1), ("G1\td__A\nG2\td__A\textra\n", 2)] {
        assert!(
            matches!(read(contents), Err(CuratomaticError::TaxonomyFile { line: l, .. }) if l == line),
            "{contents:?}"
        );
    }
}

#[test]
fn an_invalid_taxonomy_is_an_error() {
    assert!(matches!(
        read("G1\td__A;x__B\n"),
        Err(CuratomaticError::TaxonomyFile { line: 1, .. })
    ));
}

#[test]
fn a_genome_given_twice_is_an_error() {
    assert!(matches!(
        read("G000000001\td__A\nGB_GCA_000000001.1\td__A\n"),
        Err(CuratomaticError::TaxonomyFile { line: 2, .. })
    ));
}

#[test]
fn compare_finds_the_missing_and_extra_genomes() {
    let file = read("G1\td__A\nG2\td__A\n").unwrap();
    let gids = BTreeSet::from([GenomeId::new("G2"), GenomeId::new("G3")]);
    assert_eq!(
        file.compare(&gids),
        (
            BTreeSet::from([GenomeId::new("G3")]),
            BTreeSet::from([GenomeId::new("G1")])
        )
    );
}

#[test]
fn written_rows_match_the_input_format_and_read_back() {
    let gid = GenomeId::new("G003697165");
    let taxa: Vec<Taxon> = ECOLI.split(';').map(Taxon::new).collect();
    let tsv = to_tsv([(&gid, taxa.as_slice())]);
    assert_eq!(tsv, format!("G003697165\t{ECOLI}\n"));

    let file = read(&tsv).unwrap();
    assert_eq!(file.get(&gid), Some(&Taxonomy::from_string(ECOLI).unwrap()));
}

#[test]
fn from_file_reads_a_missing_file_as_an_io_error() {
    assert!(matches!(
        TaxonomyFile::from_file(Path::new("does/not/exist.tsv")),
        Err(CuratomaticError::Io(_))
    ));
}
