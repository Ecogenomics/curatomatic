use std::collections::BTreeSet;
use std::fmt;

use strum::IntoEnumIterator;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::rank::Rank;
use crate::taxon::Taxon;

/// The taxa assigned to each rank.
///
/// A rank may hold more than one taxon (e.g. when a name is shared), or none at all. Only named
/// taxa are held: a placeholder (e.g. `g__`) is refused by [`Taxonomy::add_taxon`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Taxonomy {
    pub(crate) d: Vec<Taxon>,
    pub(crate) p: Vec<Taxon>,
    pub(crate) c: Vec<Taxon>,
    pub(crate) o: Vec<Taxon>,
    pub(crate) f: Vec<Taxon>,
    pub(crate) g: Vec<Taxon>,
    pub(crate) s: Vec<Taxon>,
}

impl Taxonomy {
    /// Parse a semicolon separated taxonomy string (e.g. `d__Bacteria; p__Pseudomonadota`). A
    /// placeholder is an error, see [`Taxonomy::add_taxon`].
    pub(crate) fn from_string(tax_string: &str) -> CuratomaticResult<Self> {
        let mut out = Taxonomy::default();
        for taxon in tax_string.split(';') {
            out.add_taxon(Taxon::new(taxon.trim()))?;
        }
        Ok(out)
    }

    /// Iterate over every taxon, ordered from the highest rank to the lowest.
    pub(crate) fn iter_ranks(&self) -> impl Iterator<Item = &Taxon> {
        Rank::iter().flat_map(|rank| self.get_rank(rank).iter())
    }

    /// Remove a taxon from the rank indicated by its prefix.
    pub(crate) fn remove_taxon(&mut self, taxon: &Taxon) -> CuratomaticResult<()> {
        let rank = taxon.rank()?;
        let taxa = self.rank_mut(rank);
        let idx = taxa
            .iter()
            .position(|x| x == taxon)
            .ok_or_else(|| CuratomaticError::TaxonNotFound(taxon.to_string()))?;
        taxa.remove(idx);
        Ok(())
    }

    /// Add a taxon to the rank indicated by its prefix. A malformed taxon is an error (see
    /// [`Taxon::rank`]), and so is a placeholder ([`CuratomaticError::Placeholder`]).
    pub(crate) fn add_taxon(&mut self, taxon: Taxon) -> CuratomaticResult<()> {
        let rank = taxon.rank()?;
        if !taxon.is_named() {
            return Err(CuratomaticError::Placeholder(taxon.to_string()));
        }
        self.rank_mut(rank).push(taxon);
        Ok(())
    }

    /// True if every rank has at least one taxon assigned.
    pub(crate) fn all_assigned(&self) -> bool {
        Rank::iter().all(|rank| !self.get_rank(rank).is_empty())
    }

    /// The ranks that have no taxon assigned.
    pub(crate) fn missing_ranks(&self) -> BTreeSet<Rank> {
        Rank::iter()
            .filter(|&rank| self.get_rank(rank).is_empty())
            .collect()
    }

    /// The taxa assigned to a rank.
    pub(crate) fn get_rank(&self, rank: Rank) -> &[Taxon] {
        match rank {
            Rank::Domain => &self.d,
            Rank::Phylum => &self.p,
            Rank::Class => &self.c,
            Rank::Order => &self.o,
            Rank::Family => &self.f,
            Rank::Genus => &self.g,
            Rank::Species => &self.s,
        }
    }

    fn rank_mut(&mut self, rank: Rank) -> &mut Vec<Taxon> {
        match rank {
            Rank::Domain => &mut self.d,
            Rank::Phylum => &mut self.p,
            Rank::Class => &mut self.c,
            Rank::Order => &mut self.o,
            Rank::Family => &mut self.f,
            Rank::Genus => &mut self.g,
            Rank::Species => &mut self.s,
        }
    }
}

/// Formats as a `; ` separated string, ordered from the highest rank to the lowest.
impl fmt::Display for Taxonomy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, taxon) in self.iter_ranks().enumerate() {
            if i > 0 {
                f.write_str("; ")?;
            }
            write!(f, "{taxon}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
