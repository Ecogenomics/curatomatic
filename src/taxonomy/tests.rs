use super::*;

const FULL: &str = "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; \
                    o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; \
                    s__Escherichia coli";

#[test]
fn round_trip() {
    let tax = Taxonomy::from_string(FULL).unwrap();
    assert_eq!(tax.to_string(), FULL);
    assert!(tax.all_assigned());
    assert!(tax.missing_ranks().is_empty());
}

#[test]
fn parses_without_spaces_and_out_of_order() {
    let tax = Taxonomy::from_string("g__Escherichia;d__Bacteria").unwrap();
    assert_eq!(tax.to_string(), "d__Bacteria; g__Escherichia");
    assert_eq!(tax.get_rank(Rank::Genus), [Taxon::new("g__Escherichia")]);
}

#[test]
fn missing_ranks() {
    let tax = Taxonomy::from_string("d__Bacteria; p__Pseudomonadota").unwrap();
    assert!(!tax.all_assigned());
    assert_eq!(
        tax.missing_ranks().into_iter().collect::<Vec<_>>(),
        vec![
            Rank::Class,
            Rank::Order,
            Rank::Family,
            Rank::Genus,
            Rank::Species
        ]
    );
}

#[test]
fn add_and_remove() {
    let mut tax = Taxonomy::default();
    tax.add_taxon(Taxon::new("g__Foo")).unwrap();
    tax.add_taxon(Taxon::new("g__Bar")).unwrap();
    assert_eq!(
        tax.iter_ranks().map(Taxon::as_str).collect::<Vec<_>>(),
        ["g__Foo", "g__Bar"]
    );
    tax.remove_taxon(&Taxon::new("g__Foo")).unwrap();
    assert_eq!(
        tax.iter_ranks().map(Taxon::as_str).collect::<Vec<_>>(),
        ["g__Bar"]
    );
    assert_eq!(
        tax.remove_taxon(&Taxon::new("g__Foo")),
        Err(CuratomaticError::TaxonNotFound("g__Foo".to_owned()))
    );
}

#[test]
fn invalid_taxa() {
    assert_eq!(Taxonomy::from_string(""), Err(CuratomaticError::EmptyTaxon));
    assert_eq!(
        Taxonomy::from_string("d__Bacteria;"),
        Err(CuratomaticError::EmptyTaxon)
    );
    assert_eq!(
        Taxonomy::from_string("x__Foo"),
        Err(CuratomaticError::UnknownRank("x__Foo".to_owned()))
    );
}

#[test]
fn a_placeholder_is_refused() {
    let mut taxonomy = Taxonomy::default();
    assert_eq!(
        taxonomy.add_taxon(Taxon::new("g__")),
        Err(CuratomaticError::Placeholder("g__".to_owned()))
    );
    assert_eq!(taxonomy, Taxonomy::default());
    assert_eq!(
        Taxonomy::from_string("d__Bacteria;p__"),
        Err(CuratomaticError::Placeholder("p__".to_owned()))
    );
}

// Expected values for the duplicate rank and ordering tests were produced by the Python
// implementation (`curatomatic/model/taxonomy.py`, `Taxonomy.from_string` and `iter_ranks`):
// `Taxonomy.from_string('p__B; g__Z; p__A; d__X; g__Y')` has `to_string()`
// `d__X; p__B; p__A; g__Z; g__Y`.
#[test]
fn duplicate_ranks_keep_every_taxon_in_input_order() {
    let tax = Taxonomy::from_string("p__A; p__B").unwrap();
    assert_eq!(
        tax.get_rank(Rank::Phylum),
        [Taxon::new("p__A"), Taxon::new("p__B")]
    );
    assert_eq!(tax.to_string(), "p__A; p__B");
}

#[test]
fn iter_ranks_orders_by_rank_then_by_input_order() {
    let tax = Taxonomy::from_string("p__B; g__Z; p__A; d__X; g__Y").unwrap();
    assert_eq!(
        tax.iter_ranks().map(Taxon::as_str).collect::<Vec<_>>(),
        ["d__X", "p__B", "p__A", "g__Z", "g__Y"]
    );
    assert_eq!(tax.to_string(), "d__X; p__B; p__A; g__Z; g__Y");
}

#[test]
fn display_joins_with_a_semicolon_and_space() {
    let tax = Taxonomy::from_string("d__Bacteria;g__Escherichia").unwrap();
    assert_eq!(tax.to_string(), "d__Bacteria; g__Escherichia");
    assert_eq!(Taxonomy::default().to_string(), "");
}

#[test]
fn empty_taxonomy_has_every_rank_missing() {
    let tax = Taxonomy::default();
    assert!(!tax.all_assigned());
    assert_eq!(tax.missing_ranks().len(), Rank::iter().count());
    assert_eq!(tax.iter_ranks().count(), 0);
}

#[test]
fn remove_taxon_removes_only_the_first_match_in_its_rank() {
    let mut tax = Taxonomy::from_string("g__A; g__B; g__A").unwrap();
    tax.remove_taxon(&Taxon::new("g__A")).unwrap();
    assert_eq!(tax.to_string(), "g__B; g__A");
}

#[test]
fn remove_taxon_with_an_unknown_rank_is_an_error() {
    let mut tax = Taxonomy::default();
    assert_eq!(
        tax.remove_taxon(&Taxon::new("x__A")),
        Err(CuratomaticError::UnknownRank("x__A".to_owned()))
    );
}

#[test]
fn spaces_around_taxa_are_trimmed_but_inner_spaces_are_kept() {
    let tax = Taxonomy::from_string("  s__Escherichia coli  ;  g__Escherichia ").unwrap();
    assert_eq!(tax.to_string(), "g__Escherichia; s__Escherichia coli");
}

#[test]
fn a_taxonomy_with_a_taxon_without_a_double_underscore_is_an_error() {
    assert_eq!(
        Taxonomy::from_string("d__Bacteria;pFoo"),
        Err(CuratomaticError::InvalidTaxon("pFoo".to_owned()))
    );
}
