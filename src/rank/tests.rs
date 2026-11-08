use super::*;

#[test]
fn prefix_round_trips_for_every_rank() {
    for rank in Rank::iter() {
        assert_eq!(Rank::try_from(rank.prefix()).ok(), Some(rank));
    }
}

#[test]
fn iter_is_strictly_ordered_from_highest_to_lowest() {
    let ranks: Vec<Rank> = Rank::iter().collect();
    assert!(ranks.windows(2).all(|pair| pair[0] < pair[1]));
}

#[test]
fn unknown_prefix_is_rejected() {
    assert!(matches!(
        Rank::try_from('x'),
        Err(CuratomaticError::UnknownRank(found)) if found == "x"
    ));
}

#[test]
fn of_taxon_reads_the_leading_prefix() {
    assert_eq!(Rank::of_taxon("g__Escherichia").ok(), Some(Rank::Genus));
    assert!(matches!(
        Rank::of_taxon("x__Foo"),
        Err(CuratomaticError::UnknownRank(found)) if found == "x__Foo"
    ));
    assert!(matches!(
        Rank::of_taxon(""),
        Err(CuratomaticError::EmptyTaxon)
    ));
}

#[test]
fn parent_of_the_domain_wraps_around_to_the_species() {
    // Python looks the parent up with `RANKS[RANKS.index(rank) - 1]`, and `RANKS[-1]` is the
    // species (`curatomatic/model/rank.py`).
    assert_eq!(Rank::Domain.parent(), Rank::Species);
}

#[test]
fn parent_is_the_next_higher_rank() {
    assert_eq!(Rank::Phylum.parent(), Rank::Domain);
    assert_eq!(Rank::Class.parent(), Rank::Phylum);
    assert_eq!(Rank::Order.parent(), Rank::Class);
    assert_eq!(Rank::Family.parent(), Rank::Order);
    assert_eq!(Rank::Genus.parent(), Rank::Family);
    assert_eq!(Rank::Species.parent(), Rank::Genus);
}

#[test]
fn child_is_the_next_lower_rank_and_the_species_has_none() {
    assert_eq!(Rank::Domain.child(), Some(Rank::Phylum));
    assert_eq!(Rank::Phylum.child(), Some(Rank::Class));
    assert_eq!(Rank::Class.child(), Some(Rank::Order));
    assert_eq!(Rank::Order.child(), Some(Rank::Family));
    assert_eq!(Rank::Family.child(), Some(Rank::Genus));
    assert_eq!(Rank::Genus.child(), Some(Rank::Species));
    assert_eq!(Rank::Species.child(), None);
}

#[test]
fn prefixes_are_the_gtdb_rank_letters() {
    let prefixes: String = Rank::iter().map(Rank::prefix).collect();
    assert_eq!(prefixes, "dpcofgs");
}

#[test]
fn try_from_is_case_sensitive() {
    assert!(Rank::try_from('G').is_err());
}
