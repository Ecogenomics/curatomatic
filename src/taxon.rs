use std::collections::BTreeSet;
use std::fmt;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::metadata::Metadata;
use crate::rank::Rank;

/// Characters that are easily confused with one another when a name is read by a person.
const BAD_PAIRS: [(char, char); 5] = [('I', 'l'), ('0', 'O'), ('1', 'I'), ('5', 'S'), ('2', 'Z')];

/// The name of a single taxon, including its rank prefix (e.g. `g__Escherichia`).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Taxon(String);

impl Taxon {
    pub(crate) fn new(name: impl Into<String>) -> Self {
        Taxon(name.into())
    }

    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }

    /// The rank of this taxon, determined by the first character of its name.
    pub(crate) fn rank(&self) -> CuratomaticResult<Rank> {
        Rank::of_taxon(self.as_str())
    }
}

impl fmt::Display for Taxon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl AsRef<str> for Taxon {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl From<String> for Taxon {
    fn from(name: String) -> Self {
        Taxon(name)
    }
}

impl From<&str> for Taxon {
    fn from(name: &str) -> Self {
        Taxon(name.to_owned())
    }
}

/// `1` if the taxon contains a pair of easily confused characters, otherwise `0`.
fn score_candidate_taxon(taxon: &Taxon) -> u8 {
    let name = taxon.as_str();
    u8::from(
        BAD_PAIRS
            .iter()
            .any(|&(a, b)| name.contains(a) && name.contains(b)),
    )
}

/// Generate a taxon for `rank` that is not already in `known_taxa`.
///
/// The NCBI WGS prefix of the genomes is used when there is one, preferring the most readable. If
/// none is available then the name is derived from the lowest genome ID.
pub(crate) fn generate_unique_rank(
    rank: Rank,
    leaf_node_gids: &BTreeSet<GenomeId>,
    metadata: &Metadata,
    known_taxa: &BTreeSet<Taxon>,
) -> CuratomaticResult<Taxon> {
    let wgs_available: BTreeSet<Taxon> = leaf_node_gids
        .iter()
        .filter_map(|gid| metadata.rows.get(gid)?.wgs.as_ref())
        .map(|wgs| Taxon::new(format!("{}__{wgs}", rank.prefix())))
        .filter(|taxon| !known_taxa.contains(taxon))
        .collect();

    // Score the WGS sequences by how readable they are, the sets are sorted so the first of the
    // best scoring is also the lowest.
    if let Some(taxon) = wgs_available.iter().min_by_key(|taxon| score_candidate_taxon(taxon)) {
        return Ok(taxon.clone());
    }

    let gid = leaf_node_gids
        .first()
        .ok_or_else(|| CuratomaticError::Invariant("no genomes to name a taxon after".to_owned()))?;
    Ok(Taxon::new(format!("{}__NOVEL-{gid}", rank.prefix())))
}

#[cfg(test)]
mod tests;
