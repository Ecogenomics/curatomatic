use std::fmt;

/// Leaf nodes with this prefix are zombies: copies of a genome that is already in the tree.
const ZOMBIE_PREFIX: &str = "D-";

/// The ID of a genome, e.g. `G002502135`.
///
/// The leaf of a zombie (e.g. `D-G002502135`) also holds a `GenomeId`, see
/// [`GenomeId::is_zombie`]. The derived order is byte-wise, which is the same as the plain string
/// order that Python sorts genome IDs by (`util/taxon.py:26`).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct GenomeId(String);

impl GenomeId {
    pub(crate) fn new(gid: impl Into<String>) -> Self {
        GenomeId(gid.into())
    }

    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }

    /// True if this is the ID of a zombie leaf: a copy of a genome that is already in the tree.
    pub(crate) fn is_zombie(&self) -> bool {
        self.0.starts_with(ZOMBIE_PREFIX)
    }
}

impl fmt::Display for GenomeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Convert a genome ID to canonical format, e.g. `GCA_947480885.1` -> `G947480885`.
///
/// The prefix (`GB_GCA_`, `RS_GCF_`, `GCA_` or `GCF_`) and the two character version suffix are
/// removed. IDs with any other prefix are returned unchanged.
pub(crate) fn gid_to_canonical(gid: &str) -> GenomeId {
    let Some(rest) = ["GB_GCA_", "RS_GCF_", "GCA_", "GCF_"]
        .into_iter()
        .find_map(|prefix| gid.strip_prefix(prefix))
    else {
        return GenomeId::new(gid);
    };

    // Drop the last two characters (the `.1` version suffix).
    let end = rest.char_indices().rev().nth(1).map_or(0, |(idx, _)| idx);
    GenomeId::new(format!("G{}", rest.get(..end).unwrap_or_default()))
}

#[cfg(test)]
mod tests;
