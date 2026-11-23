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
fn parent_is_the_next_higher_rank_and_the_domain_has_none() {
    // The domain must not wrap around to the species.
    assert_eq!(Rank::Domain.parent(), None);
    assert_eq!(Rank::Phylum.parent(), Some(Rank::Domain));
    assert_eq!(Rank::Class.parent(), Some(Rank::Phylum));
    assert_eq!(Rank::Order.parent(), Some(Rank::Class));
    assert_eq!(Rank::Family.parent(), Some(Rank::Order));
    assert_eq!(Rank::Genus.parent(), Some(Rank::Family));
    assert_eq!(Rank::Species.parent(), Some(Rank::Genus));
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

#[test]
fn from_name_ignores_case_and_round_trips() {
    for rank in Rank::iter() {
        assert_eq!(Rank::from_name(rank.name()), Ok(rank));
        assert_eq!(Rank::from_name(&rank.name().to_lowercase()), Ok(rank));
        assert_eq!(Rank::from_name(&rank.name().to_uppercase()), Ok(rank));
    }
    assert_eq!(Rank::from_name(" Genus "), Ok(Rank::Genus));
}

#[test]
fn from_name_rejects_prefixes_and_other_names() {
    for name in ["p", "Phyla", "Phylumx", "Kingdom", ""] {
        assert_eq!(
            Rank::from_name(name),
            Err(CuratomaticError::UnknownRank(name.to_owned()))
        );
    }
}

#[test]
fn a_taxon_without_a_double_underscore_after_its_prefix_is_invalid() {
    for taxon in ["gFoo", "g_Foo", "g", "g_"] {
        assert_eq!(
            Rank::of_taxon(taxon),
            Err(CuratomaticError::InvalidTaxon(taxon.to_owned())),
            "{taxon}"
        );
    }
}

#[test]
fn a_name_starting_with_an_underscore_or_space_is_invalid() {
    // A trailing space is invalid too: a name is never trimmed or otherwise edited.
    for taxon in [
        "g___Foo",
        "g__ Foo",
        "g__Foo ",
        "g__\tFoo",
        "s__Escherichia coli\t",
    ] {
        assert_eq!(
            Rank::of_taxon(taxon),
            Err(CuratomaticError::InvalidTaxon(taxon.to_owned())),
            "{taxon:?}"
        );
    }
}

#[test]
fn every_prefix_followed_by_a_double_underscore_has_its_rank() {
    for rank in Rank::iter() {
        let prefix = rank.prefix();
        for name in [
            "Name",
            "Name_A",
            "CAG-698",
            "Name two",
            "{unresolved} sp1",
            "",
        ] {
            let taxon = format!("{prefix}__{name}");
            assert_eq!(Rank::of_taxon(&taxon), Ok(rank), "{taxon}");
            assert_eq!(Rank::split_taxon(&taxon), Ok((rank, name)), "{taxon}");
        }
    }
}

#[test]
fn an_unknown_prefix_with_a_double_underscore_is_an_unknown_rank() {
    for taxon in ["x__Foo", "G__Foo"] {
        assert_eq!(
            Rank::of_taxon(taxon),
            Err(CuratomaticError::UnknownRank(taxon.to_owned())),
            "{taxon}"
        );
    }
}
