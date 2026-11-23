use strum::{EnumIter, IntoEnumIterator};

use crate::error::{CuratomaticError, CuratomaticResult};

/// Separates the rank prefix of a taxon from its name (e.g. `g__Escherichia`).
pub(crate) const PREFIX_SEPARATOR: &str = "__";

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

    /// The full name of this rank, capitalised (e.g. `Phylum`).
    pub(crate) fn name(self) -> &'static str {
        match self {
            Rank::Domain => "Domain",
            Rank::Phylum => "Phylum",
            Rank::Class => "Class",
            Rank::Order => "Order",
            Rank::Family => "Family",
            Rank::Genus => "Genus",
            Rank::Species => "Species",
        }
    }

    /// The rank with this full name, ignoring case (e.g. `Phylum` or `phylum`).
    pub(crate) fn from_name(name: &str) -> CuratomaticResult<Rank> {
        let name = name.trim();
        Rank::iter()
            .find(|rank| rank.name().eq_ignore_ascii_case(name))
            .ok_or_else(|| CuratomaticError::UnknownRank(name.to_owned()))
    }

    /// The next higher rank, or `None` for the domain.
    pub(crate) fn parent(self) -> Option<Rank> {
        Rank::iter().rev().skip_while(|&rank| rank != self).nth(1)
    }

    /// The next lower rank, or `None` for the species.
    pub(crate) fn child(self) -> Option<Rank> {
        Rank::iter().skip_while(|&rank| rank != self).nth(1)
    }

    /// The rank of a taxon, determined by its prefix (e.g. `g__Escherichia`). See
    /// [`Rank::split_taxon`] for the format of a taxon.
    pub(crate) fn of_taxon(taxon: &str) -> CuratomaticResult<Rank> {
        Rank::split_taxon(taxon).map(|(rank, _name)| rank)
    }

    /// Split a taxon into its rank and its name: `g__Escherichia` is the genus `Escherichia`, and
    /// a placeholder (`g__`) is a genus with an empty name.
    ///
    /// A taxon is a rank prefix, then [`PREFIX_SEPARATOR`], then its name. The name must not start
    /// with `_` or whitespace, and must not end with whitespace. Anything else in the name is
    /// allowed, and is never edited.
    ///
    /// An empty taxon is [`CuratomaticError::EmptyTaxon`]. A well formed taxon with a prefix that
    /// is not a rank (`x__Foo`, `G__Foo`) is [`CuratomaticError::UnknownRank`]. Any other
    /// malformed taxon (`gFoo`, `g_Foo`, `g___Foo`, `g__ Foo`, `g__Foo `) is
    /// [`CuratomaticError::InvalidTaxon`].
    pub(crate) fn split_taxon(taxon: &str) -> CuratomaticResult<(Rank, &str)> {
        let invalid = || CuratomaticError::InvalidTaxon(taxon.to_owned());
        let mut chars = taxon.chars();
        let prefix = chars.next().ok_or(CuratomaticError::EmptyTaxon)?;
        let name = chars
            .as_str()
            .strip_prefix(PREFIX_SEPARATOR)
            .ok_or_else(invalid)?;
        let rank = Rank::try_from(prefix)
            .map_err(|_unknown_prefix| CuratomaticError::UnknownRank(taxon.to_owned()))?;
        let bad_start = name
            .chars()
            .next()
            .is_some_and(|c| c == '_' || c.is_whitespace());
        let bad_end = name.chars().next_back().is_some_and(char::is_whitespace);
        if bad_start || bad_end {
            return Err(invalid());
        }
        Ok((rank, name))
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
