use strum::{EnumIter, IntoEnumIterator};

use crate::error::{CuratomaticError, CuratomaticResult};

/// A taxonomic rank, in order from highest (domain) to lowest (species).
///
/// Iterate over every rank, highest to lowest, with [`Rank::iter`] (from [`strum::IntoEnumIterator`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, EnumIter)]
pub(crate) enum Rank {
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
}

impl Rank {
    /// The single character prefix used to denote this rank (e.g. `d` for `d__Bacteria`).
    pub(crate) fn prefix(self) -> char {
        // Exhaustive, with no catch-all: adding a variant fails to compile until it is given a
        // prefix here. Parsing (`TryFrom<char>`) is derived from this, so it cannot drift.
        match self {
            Rank::Domain => 'd',
            Rank::Phylum => 'p',
            Rank::Class => 'c',
            Rank::Order => 'o',
            Rank::Family => 'f',
            Rank::Genus => 'g',
            Rank::Species => 's',
        }
    }

    /// The next higher rank, or `None` for the domain.
    pub(crate) fn parent(self) -> Option<Rank> {
        Rank::iter().rev().skip_while(|&rank| rank != self).nth(1)
    }

    /// The next lower rank, or `None` for the species.
    pub(crate) fn child(self) -> Option<Rank> {
        Rank::iter().skip_while(|&rank| rank != self).nth(1)
    }

    /// The rank of a taxon, determined by the first character of its name (e.g. `g__Escherichia`).
    pub(crate) fn of_taxon(taxon: &str) -> CuratomaticResult<Rank> {
        let prefix = taxon.chars().next().ok_or(CuratomaticError::EmptyTaxon)?;
        Rank::try_from(prefix)
            .map_err(|_unknown_prefix| CuratomaticError::UnknownRank(taxon.to_owned()))
    }
}

impl TryFrom<char> for Rank {
    type Error = CuratomaticError;

    /// The rank denoted by a prefix character, if it is a known one.
    fn try_from(prefix: char) -> Result<Self, Self::Error> {
        Rank::iter()
            .find(|rank| rank.prefix() == prefix)
            .ok_or_else(|| CuratomaticError::UnknownRank(prefix.to_string()))
    }
}

#[cfg(test)]
mod tests;
