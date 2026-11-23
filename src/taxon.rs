use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::metadata::Metadata;
use crate::rank::Rank;

/// Characters that are easily confused with one another when a name is read by a person.
const BAD_PAIRS: [(char, char); 5] = [('I', 'l'), ('0', 'O'), ('1', 'I'), ('5', 'S'), ('2', 'Z')];

/// The name of a novel taxon that could not be given a unique name, after its rank prefix.
const UNRESOLVABLE: &str = "UNRESOLVABLE";

/// The start of a species name that curatomatic must not change, e.g. `s__{unresolved} NAME`.
const UNRESOLVED_SPECIES: &str = "s__{unresolved}";

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

    /// The rank of this taxon, determined by its prefix. A malformed taxon is an error, see
    /// [`Rank::split_taxon`].
    pub(crate) fn rank(&self) -> CuratomaticResult<Rank> {
        Rank::of_taxon(self.as_str())
    }

    /// The sentinel for a novel taxon of `rank` that could not be given a unique name, e.g.
    /// `g__UNRESOLVABLE`. Unlike every other taxon, it may be on more than one node.
    pub(crate) fn unresolvable(rank: Rank) -> Self {
        Taxon(format!("{}__{UNRESOLVABLE}", rank.prefix()))
    }

    /// True if this is the sentinel of [`Taxon::unresolvable`].
    pub(crate) fn is_unresolvable(&self) -> bool {
        self.rank()
            .is_ok_and(|rank| self == &Taxon::unresolvable(rank))
    }

    /// True if this is a name: a valid taxon with text after its `<prefix>__` (see
    /// [`Rank::split_taxon`]), rather than a placeholder that is only the prefix (e.g. `g__`).
    pub(crate) fn is_named(&self) -> bool {
        Rank::split_taxon(self.as_str()).is_ok_and(|(_rank, name)| !name.is_empty())
    }

    /// True if this is a species of the form `s__{unresolved} NAME`.
    pub(crate) fn is_unresolved_species(&self) -> bool {
        self.as_str().starts_with(UNRESOLVED_SPECIES)
    }

    /// True if this is a valid name: named (see [`Taxon::is_named`]), and not an unresolved
    /// species (see [`Taxon::is_unresolved_species`]). A novel taxon is never created above one.
    pub(crate) fn is_valid_name(&self) -> bool {
        self.is_named() && !self.is_unresolved_species()
    }

    /// True if this is a Latin name: a capital letter followed by lower case letters, with an
    /// optional suffix of an underscore and capital letters (`g__Escherichia`, `g__Luteimonas_D`).
    /// Placeholder style names (`g__JANWYU01`, `g__2-01-FULL-41-22`, `g__CAG-698`) are not.
    pub(crate) fn is_latin_name(&self) -> bool {
        let Some((_, name)) = self.as_str().split_once("__") else {
            return false;
        };
        let (stem, suffix) = name.split_once('_').unwrap_or((name, ""));
        let mut chars = stem.chars();
        chars.next().is_some_and(|c| c.is_ascii_uppercase())
            && !chars.as_str().is_empty()
            && chars.all(|c| c.is_ascii_lowercase())
            && (!name.contains('_')
                || (!suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_uppercase())))
    }

    /// The genus that a binomial species name puts its genome in: `g__Escherichia` for
    /// `s__Escherichia coli` (or `s__Escherichia sp000123`). `None` for anything that is not a
    /// binomial species, including unresolved species (`s__{unresolved} NAME`).
    pub(crate) fn genus_of_species(&self) -> Option<Taxon> {
        if self.is_unresolved_species() {
            return None;
        }
        let (genus, epithet) = self.as_str().strip_prefix("s__")?.split_once(' ')?;
        (!genus.is_empty() && !epithet.trim().is_empty()).then(|| Taxon(format!("g__{genus}")))
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

/// The number of nodes that hold each named taxon (placeholders, e.g. `g__`, are not counted).
///
/// A taxon is on one node only, except [`Taxon::unresolvable`], so a count is above `1` only for
/// that sentinel, or when a taxon has been copied instead of moved.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct TaxonCounts(BTreeMap<Taxon, usize>);

impl TaxonCounts {
    /// Count each named taxon, one for each time it is given.
    pub(crate) fn from_taxa<'a>(taxa: impl IntoIterator<Item = &'a Taxon>) -> Self {
        let mut out = BTreeMap::new();
        for taxon in taxa.into_iter().filter(|taxon| taxon.is_named()) {
            let n: &mut usize = out.entry(taxon.clone()).or_default();
            *n = n.saturating_add(1);
        }
        TaxonCounts(out)
    }

    /// The number of nodes that hold `taxon`, `0` if none do.
    pub(crate) fn get(&self, taxon: &Taxon) -> usize {
        self.0.get(taxon).copied().unwrap_or_default()
    }

    /// Every taxon and its count, ordered by taxon.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&Taxon, usize)> {
        self.0.iter().map(|(taxon, &n)| (taxon, n))
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
/// none is available then the taxon is named after the lowest genome ID whose name is not already
/// known (e.g. `g__G963674675`). If every genome ID is taken (e.g. the input tree named a taxon
/// after a genome that is not in it), the taxon is the sentinel [`Taxon::unresolvable`].
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
    if let Some(taxon) = wgs_available
        .iter()
        .min_by_key(|taxon| score_candidate_taxon(taxon))
    {
        return Ok(taxon.clone());
    }

    // Genome IDs are sorted, so this tries the lowest first.
    let gid_taxon = |gid: &GenomeId| Taxon::new(format!("{}__{gid}", rank.prefix()));
    if let Some(taxon) = leaf_node_gids
        .iter()
        .map(gid_taxon)
        .find(|taxon| !known_taxa.contains(taxon))
    {
        return Ok(taxon);
    }

    if leaf_node_gids.is_empty() {
        return Err(CuratomaticError::Invariant(
            "no genomes to name a taxon after".to_owned(),
        ));
    }
    Ok(Taxon::unresolvable(rank))
}

#[cfg(test)]
mod tests;
