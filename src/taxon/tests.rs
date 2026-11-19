use super::*;
use std::fmt::Write as _;
use std::io::Cursor;
use std::path::Path;

// Expected values for `score_candidate_taxon` and the WGS choice in `generate_unique_rank` come from
// the original Python implementation, run on the same inputs. The genome ID fallback differs on
// purpose: curators name these taxa after the bare genome ID, not `NOVEL-<gid>`.

fn gids(gids: &[&str]) -> BTreeSet<GenomeId> {
    gids.iter().map(|&gid| GenomeId::new(gid)).collect()
}

fn known(taxa: &[&str]) -> BTreeSet<Taxon> {
    taxa.iter().map(|&taxon| Taxon::new(taxon)).collect()
}

/// Metadata where each genome is given as `(accession, wgs)`.
fn metadata(rows: &[(&str, &str)]) -> Metadata {
    let mut tsv = String::from("accession\tncbi_wgs_formatted\n");
    for (accession, wgs) in rows {
        writeln!(tsv, "{accession}\t{wgs}").unwrap();
    }
    Metadata::from_reader(
        Cursor::new(tsv),
        Path::new("names.tsv"),
        "ncbi_wgs_formatted",
    )
    .unwrap()
}

fn generate(
    rank: Rank,
    leaf_gids: &[&str],
    metadata: &Metadata,
    known_taxa: &[&str],
) -> CuratomaticResult<Taxon> {
    generate_unique_rank(rank, &gids(leaf_gids), metadata, &known(known_taxa))
}

#[test]
fn score_is_one_for_each_bad_pair_in_either_order() {
    for (a, b) in BAD_PAIRS {
        for name in [
            format!("g__{a}{b}"),
            format!("g__{b}{a}"),
            format!("g__x{a}y{b}z"),
        ] {
            assert_eq!(
                score_candidate_taxon(&Taxon::new(name.clone())),
                1,
                "{name}"
            );
        }
    }
}

#[test]
fn score_is_zero_without_a_complete_bad_pair() {
    for name in ["g__abc", "g__I", "g__1", "g__SS", "g__l", "g__O", "g__"] {
        assert_eq!(score_candidate_taxon(&Taxon::new(name)), 0, "{name}");
    }
}

#[test]
fn score_is_one_when_any_pair_is_present() {
    // `I`+`l` and `0`+`O` are both present, only one pair is needed.
    assert_eq!(score_candidate_taxon(&Taxon::new("g__Il0")), 1);
    assert_eq!(score_candidate_taxon(&Taxon::new("g__Z2x")), 1);
}

#[test]
fn without_metadata_the_lowest_gid_names_the_taxon() {
    let taxon = generate(Rank::Genus, &["G3", "G1", "G2"], &Metadata::default(), &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__G1"));
}

#[test]
fn the_lowest_gid_is_chosen_by_plain_string_order() {
    // Python `sorted(...)[0]`: `G10` sorts before `G9`.
    let taxon = generate(Rank::Genus, &["G9", "G10"], &Metadata::default(), &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__G10"));
}

#[test]
fn the_rank_prefix_follows_the_rank() {
    let taxon = generate(Rank::Species, &["G1"], &Metadata::default(), &[]).unwrap();
    assert_eq!(taxon, Taxon::new("s__G1"));
    let meta = metadata(&[("G1", "ABC")]);
    assert_eq!(
        generate(Rank::Family, &["G1"], &meta, &[]).unwrap(),
        Taxon::new("f__ABC")
    );
}

#[test]
fn wgs_is_used_when_available() {
    let meta = metadata(&[("G1", "ABC")]);
    let taxon = generate(Rank::Genus, &["G1"], &meta, &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__ABC"));
}

#[test]
fn a_wgs_that_is_already_known_is_skipped() {
    let meta = metadata(&[("G1", "ABC"), ("G2", "DEF")]);
    let taxon = generate(Rank::Genus, &["G1", "G2"], &meta, &["g__ABC"]).unwrap();
    assert_eq!(taxon, Taxon::new("g__DEF"));
}

#[test]
fn a_known_taxon_of_another_rank_does_not_conflict() {
    let meta = metadata(&[("G1", "ABC")]);
    let taxon = generate(Rank::Species, &["G1"], &meta, &["g__ABC"]).unwrap();
    assert_eq!(taxon, Taxon::new("s__ABC"));
}

#[test]
fn ties_in_score_go_to_the_lowest_string() {
    let meta = metadata(&[("G1", "ZZZ"), ("G2", "AAA"), ("G3", "MMM")]);
    let taxon = generate(Rank::Genus, &["G1", "G2", "G3"], &meta, &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__AAA"));
}

#[test]
fn the_lowest_score_wins_over_the_lowest_string() {
    // `g__AIl` has a bad pair (score 1) so `g__BBB` (score 0) wins even though it sorts later.
    let meta = metadata(&[("G1", "AIl"), ("G2", "BBB"), ("G3", "CCC")]);
    let taxon = generate(Rank::Genus, &["G1", "G2", "G3"], &meta, &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__BBB"));
}

#[test]
fn when_every_wgs_scores_the_same_the_lowest_string_wins() {
    // Both score 1, and `0` sorts before `A`.
    let meta = metadata(&[("G1", "AIl"), ("G2", "0O")]);
    let taxon = generate(Rank::Genus, &["G1", "G2"], &meta, &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__0O"));
}

#[test]
fn when_every_wgs_is_known_it_falls_back_to_the_gid() {
    let meta = metadata(&[("G1", "ABC"), ("G2", "DEF")]);
    let taxon = generate(Rank::Genus, &["G2", "G1"], &meta, &["g__ABC", "g__DEF"]).unwrap();
    assert_eq!(taxon, Taxon::new("g__G1"));
}

#[test]
fn a_known_gid_name_falls_back_to_the_next_lowest_gid() {
    let taxon = generate(
        Rank::Genus,
        &["G3", "G1", "G2"],
        &Metadata::default(),
        &["g__G1"],
    )
    .unwrap();
    assert_eq!(taxon, Taxon::new("g__G2"));
}

#[test]
fn a_gid_name_of_another_rank_does_not_conflict() {
    let taxon = generate(Rank::Genus, &["G1"], &Metadata::default(), &["f__G1"]).unwrap();
    assert_eq!(taxon, Taxon::new("g__G1"));
}

#[test]
fn when_every_gid_name_is_known_the_taxon_is_unresolvable() {
    let known_taxa = ["g__G1", "g__G2"];
    let taxon = generate(
        Rank::Genus,
        &["G2", "G1"],
        &Metadata::default(),
        &known_taxa,
    )
    .unwrap();
    assert_eq!(taxon, Taxon::new("g__UNRESOLVABLE"));
    assert!(taxon.is_unresolvable());
}

#[test]
fn only_the_sentinel_is_unresolvable() {
    for name in [
        "g__G1",
        "g__",
        "s__{unresolved} sp1",
        "f__UNRESOLVABLE_A",
        "x__UNRESOLVABLE",
    ] {
        assert!(!Taxon::new(name).is_unresolvable(), "{name}");
    }
    assert!(Taxon::new("s__UNRESOLVABLE").is_unresolvable());
}

#[test]
fn a_known_wgs_and_gid_falls_back_to_the_next_gid() {
    let meta = metadata(&[("G1", "ABC")]);
    let taxon = generate(Rank::Genus, &["G1", "G2"], &meta, &["g__ABC", "g__G1"]).unwrap();
    assert_eq!(taxon, Taxon::new("g__G2"));
}

#[test]
fn genomes_without_a_wgs_or_metadata_are_ignored() {
    // `G1` has the `none` placeholder, and `G3` is not in the metadata.
    let meta = metadata(&[("G1", "none"), ("G2", "XYZ")]);
    let taxon = generate(Rank::Genus, &["G1", "G3"], &meta, &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__G1"));
}

#[test]
fn only_the_genomes_given_are_considered() {
    let meta = metadata(&[("G1", "ABC"), ("G2", "DEF")]);
    let taxon = generate(Rank::Genus, &["G2"], &meta, &[]).unwrap();
    assert_eq!(taxon, Taxon::new("g__DEF"));
}

#[test]
fn no_genomes_is_an_error() {
    // Python raises `IndexError` from `sorted(...)[0]`.
    let meta = metadata(&[("G1", "ABC")]);
    assert!(matches!(
        generate(Rank::Genus, &[], &meta, &[]),
        Err(CuratomaticError::Invariant(_))
    ));
}

#[test]
fn taxon_helpers() {
    let taxon = Taxon::new("g__Escherichia");
    assert_eq!(taxon.as_str(), "g__Escherichia");
    assert_eq!(taxon.to_string(), "g__Escherichia");
    assert_eq!(taxon.as_ref(), "g__Escherichia");
    assert_eq!(taxon.rank(), Ok(Rank::Genus));
    assert_eq!(Taxon::from("g__Escherichia"), taxon);
    assert_eq!(Taxon::from("g__Escherichia".to_owned()), taxon);
    assert_eq!(
        Taxon::new("x__Foo").rank(),
        Err(CuratomaticError::UnknownRank("x__Foo".to_owned()))
    );
    assert_eq!(Taxon::new("").rank(), Err(CuratomaticError::EmptyTaxon));
}

#[test]
fn taxa_are_ordered_as_plain_strings() {
    assert!(Taxon::new("g__A") < Taxon::new("g__B"));
    assert!(Taxon::new("g__Z") < Taxon::new("g__a"));
}
