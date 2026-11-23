use super::*;
use strum::IntoEnumIterator;

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
fn placeholders_are_stripped_from_a_row() {
    let file = read("G1\td__B;p__;c__;o__;f__;g__;s__X y\n").unwrap();
    assert_eq!(
        file.get(&GenomeId::new("G1")).unwrap(),
        &Taxonomy::from_string("d__B;s__X y").unwrap()
    );
}

#[test]
fn a_row_of_placeholders_only_is_empty() {
    // A missing domain or species is caught later, by the checks on the leaves.
    let file = read("G1\td__;p__;g__;g__;s__\n").unwrap();
    assert_eq!(
        file.get(&GenomeId::new("G1")).unwrap(),
        &Taxonomy::default()
    );
}

#[test]
fn unresolved_species_are_kept_as_read() {
    let file = read("G1\td__Bacteria;p__;g__;s__{unresolved} sp053939565\n").unwrap();
    let taxonomy = file.get(&GenomeId::new("G1")).unwrap();
    assert_eq!(taxonomy.s, vec![Taxon::new("s__{unresolved} sp053939565")]);
}

#[test]
fn a_rank_with_a_placeholder_and_a_name_is_an_error() {
    for row in ["d__B;g__;g__X;s__X y", "d__B;g__X;g__;s__X y"] {
        let contents = format!("G1\td__B;s__A b\nG2\t{row}\n");
        assert!(
            matches!(
                read(&contents),
                Err(CuratomaticError::TaxonomyFile { line: 2, ref message, .. })
                    if message.contains("g__X")
            ),
            "{row}"
        );
    }
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

#[test]
fn a_row_with_a_malformed_taxon_is_an_error_on_its_line() {
    let error = read("G1\td__B;gFoo;s__X y\n").unwrap_err();
    assert!(
        matches!(&error, CuratomaticError::TaxonomyFile { line: 1, message, .. }
            if message.contains("gFoo")),
        "{error}"
    );
}

/// The message of the error of reading `row` as the second line, after a valid first line.
fn row_error(row: &str) -> String {
    match read(&format!("G1\td__B;s__A b\nG2\t{row}\n")) {
        Err(CuratomaticError::TaxonomyFile {
            line: 2, message, ..
        }) => message,
        other => panic!("{row}: expected an error on line 2, got {other:?}"),
    }
}

#[test]
fn a_row_out_of_rank_order_is_an_error() {
    // The file must come in rank order: it is never reordered.
    assert_eq!(
        row_error("s__A b;c__C;p__P;d__B"),
        "the taxa are not in rank order: c__C is after s__A b"
    );
    assert_eq!(
        row_error("d__B;g__G;f__F;s__A b"),
        "the taxa are not in rank order: f__F is after g__G"
    );
}

#[test]
fn a_placeholder_out_of_rank_order_is_an_error() {
    assert_eq!(
        row_error("d__B;g__;f__;s__A b"),
        "the taxa are not in rank order: f__ is after g__"
    );
}

#[test]
fn a_taxon_given_twice_in_a_row_is_an_error() {
    assert_eq!(
        row_error("d__B;o__O;o__O;s__A b"),
        "o__O is given more than once"
    );
}

#[test]
fn a_row_with_two_domains_is_an_error() {
    assert_eq!(
        row_error("d__B;d__C;s__A b"),
        "a genome has one domain only, but the row has d__B and d__C"
    );
}

#[test]
fn a_row_with_two_species_is_an_error() {
    assert_eq!(
        row_error("d__B;s__A b;s__A c"),
        "a genome has one species only, but the row has s__A b and s__A c"
    );
}

#[test]
fn several_taxa_of_a_rank_other_than_the_domain_or_species_are_kept() {
    let file = read("G1\td__B;o__O1;o__O2;f__F;g__G1;g__G2;s__A b\n").unwrap();
    let taxonomy = file.get(&GenomeId::new("G1")).unwrap();
    assert_eq!(taxonomy.o, vec![Taxon::new("o__O1"), Taxon::new("o__O2")]);
    assert_eq!(taxonomy.g, vec![Taxon::new("g__G1"), Taxon::new("g__G2")]);
}

fn read_prior(contents: &str) -> CuratomaticResult<PriorTaxonomyFile> {
    PriorTaxonomyFile::from_reader(contents.as_bytes(), Path::new("prior.tsv"))
}

/// The message of the error of reading `row` as the second line of a prior taxonomy, after a
/// valid first line.
fn prior_row_error(row: &str) -> String {
    match read_prior(&format!("G1\t{ECOLI}\nG2\t{row}\n")) {
        Err(CuratomaticError::TaxonomyFile {
            line: 2, message, ..
        }) => message,
        other => panic!("{row}: expected an error on line 2, got {other:?}"),
    }
}

/// `ECOLI` with the taxon of `rank` replaced by `taxa` (empty to drop the rank).
fn ecoli_with(rank: Rank, taxa: &[&str]) -> String {
    let mut out: Vec<String> = Vec::new();
    for taxon in ECOLI.split(';') {
        if Taxon::new(taxon).rank().unwrap() == rank {
            out.extend(taxa.iter().map(|t| (*t).to_owned()));
        } else {
            out.push(taxon.to_owned());
        }
    }
    out.join(";")
}

#[test]
fn a_full_prior_row_loads() {
    let file = read_prior(&format!(
        "G1\t{ECOLI}\nG2\td__B;p__P;c__C;o__O;f__F;g__CAG-698;s__CAG-698 sp000123\n"
    ))
    .unwrap();
    assert_eq!(
        file.get(&GenomeId::new("G1")),
        Some(&Taxonomy::from_string(ECOLI).unwrap())
    );
    assert_eq!(file.iter().count(), 2);
}

#[test]
fn a_placeholder_in_any_rank_of_a_prior_row_is_an_error() {
    for rank in Rank::iter() {
        let placeholder = format!("{}__", rank.prefix());
        let row = ecoli_with(rank, &[&placeholder]);
        assert_eq!(
            prior_row_error(&row),
            format!(
                "the {} rank is the placeholder {placeholder}, but a prior taxonomy names every \
                 rank",
                rank.name().to_lowercase()
            ),
            "{row}"
        );
    }
}

#[test]
fn a_missing_rank_in_a_prior_row_is_an_error() {
    assert_eq!(
        prior_row_error(&ecoli_with(Rank::Genus, &[])),
        "the genus rank is missing, but a prior taxonomy names every rank"
    );
    assert_eq!(
        prior_row_error("d__Bacteria"),
        "the phylum rank is missing, but a prior taxonomy names every rank"
    );
}

#[test]
fn two_taxa_of_a_rank_in_a_prior_row_is_an_error() {
    assert_eq!(
        prior_row_error(&ecoli_with(Rank::Genus, &["g__A", "g__Escherichia"])),
        "the genus rank has g__A and g__Escherichia, but a prior taxonomy has one taxon per rank"
    );
}

#[test]
fn a_prior_row_of_placeholders_only_is_an_error() {
    assert_eq!(
        prior_row_error("d__;p__;c__;o__;f__;g__;s__"),
        "the domain rank is the placeholder d__, but a prior taxonomy names every rank"
    );
}

#[test]
fn a_prior_species_that_is_not_binomial_is_an_error() {
    for species in [
        "s__{unresolved} sp1",
        "s__Escherichia",
        "s__Escherichia coli K12",
        "s__Escherichia  coli",
    ] {
        assert_eq!(
            prior_row_error(&ecoli_with(Rank::Species, &[species])),
            format!(
                "the species {species} is not a binomial name, but a prior taxonomy names every \
                 species"
            ),
        );
    }
}

#[test]
fn a_prior_species_of_another_genus_is_an_error() {
    assert_eq!(
        prior_row_error(&ecoli_with(Rank::Species, &["s__Salmonella enterica"])),
        "the species s__Salmonella enterica is not in the genus g__Escherichia"
    );
}

#[test]
fn the_shared_checks_keep_their_messages_in_a_prior_row() {
    assert_eq!(
        prior_row_error("d__B;g__G;f__F;s__G b"),
        "the taxa are not in rank order: f__F is after g__G"
    );
    assert_eq!(
        prior_row_error(&ecoli_with(
            Rank::Genus,
            &["g__Escherichia", "g__Escherichia"]
        )),
        "g__Escherichia is given more than once"
    );
}

#[test]
fn the_taxonomy_file_still_accepts_what_a_prior_rejects() {
    let rows = Rank::iter()
        .map(|rank| ecoli_with(rank, &[&format!("{}__", rank.prefix())]))
        .chain([
            ecoli_with(Rank::Genus, &[]),
            "d__Bacteria".to_owned(),
            ecoli_with(Rank::Genus, &["g__A", "g__Escherichia"]),
            "d__;p__;c__;o__;f__;g__;s__".to_owned(),
            ecoli_with(Rank::Species, &["s__{unresolved} sp1"]),
            ecoli_with(Rank::Species, &["s__Escherichia"]),
            ecoli_with(Rank::Species, &["s__Escherichia coli K12"]),
            ecoli_with(Rank::Species, &["s__Escherichia  coli"]),
            ecoli_with(Rank::Species, &["s__Salmonella enterica"]),
        ]);
    for row in rows {
        assert!(read(&format!("G1\t{row}\n")).is_ok(), "{row}");
        assert!(read_prior(&format!("G1\t{row}\n")).is_err(), "{row}");
    }
}
