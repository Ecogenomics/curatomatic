use std::collections::BTreeMap;
use std::fmt;
use std::fs;
use std::num::ParseFloatError;
use std::path::Path;
use std::str::FromStr;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::rank::Rank;

/// A relative evolutionary divergence (RED) value.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct RedValue(f64);

impl RedValue {
    /// The RED value of a leaf node.
    pub(crate) const LEAF: RedValue = RedValue(1.0);

    pub(crate) const fn new(value: f64) -> Self {
        RedValue(value)
    }

    /// The absolute difference between two RED values.
    pub(crate) const fn value(self) -> f64 {
        self.0
    }

    pub(crate) fn distance(self, other: RedValue) -> f64 {
        (self.0 - other.0).abs()
    }

    /// The RED value halfway between two RED values.
    pub(crate) fn midpoint(self, other: RedValue) -> RedValue {
        RedValue(f64::midpoint(self.0, other.0))
    }

    /// The RED value halfway from this RED value to another, as `self + (other - self) / 2`.
    ///
    /// This can differ from [`RedValue::midpoint`] in the last bit. Python places novel taxa with
    /// this form (`fill_in_missing_ranks` in `curatomatic/util/tree.py`), and a node whose RED is
    /// exactly on the threshold is only treated the same if the arithmetic is.
    pub(crate) fn halfway_to(self, other: RedValue) -> RedValue {
        RedValue(self.0 + (other.0 - self.0) / 2.0)
    }
}

impl FromStr for RedValue {
    type Err = ParseFloatError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        value.trim().parse().map(RedValue)
    }
}

/// Formats like the underlying `f64`, so a precision such as `{:.3}` is respected.
impl fmt::Display for RedValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, f)
    }
}

/// Column of the median RED file that holds the rank name.
const COL_RANK: &str = "rank";
/// Column of the median RED file that holds the median RED of the rank.
const COL_MEDIAN: &str = "median";
/// The RED value of the domain, which is the root of the tree.
const RED_DOMAIN: RedValue = RedValue::new(0.0);

/// The RED value that a rank always has, whatever the median RED file says.
fn fixed_red(rank: Rank) -> Option<RedValue> {
    match rank {
        Rank::Domain => Some(RED_DOMAIN),
        Rank::Species => Some(RedValue::LEAF),
        Rank::Phylum | Rank::Class | Rank::Order | Rank::Family | Rank::Genus => None,
    }
}

/// The median relative evolutionary divergence (RED) expected for each rank.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Red {
    data: BTreeMap<Rank, RedValue>,
}

impl Red {
    /// Read the median RED of each rank from a tab separated file, see [`Red::from_tsv`].
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        let text = fs::read_to_string(path)
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
        Self::from_tsv(&text).map_err(|e| {
            if let CuratomaticError::Red(message) = e {
                CuratomaticError::Red(format!("{}: {message}", path.display()))
            } else {
                e
            }
        })
    }

    /// Parse a tab separated file with a header and one row per rank, e.g.
    /// `median_red_rank_distribution.tsv`.
    ///
    /// The `rank` and `median` columns are found by name in the header, and any other column is
    /// ignored. A rank is its full name in any case (e.g. `Phylum`). A rank that is not given has
    /// no value, except the domain, which is always [`RED_DOMAIN`], and the species, which is
    /// always [`RedValue::LEAF`]. A row for either is accepted only if it has that value.
    ///
    /// A missing column, an unparsable median, an unknown rank and a rank given twice are errors.
    pub(crate) fn from_tsv(text: &str) -> CuratomaticResult<Self> {
        let mut lines = text
            .lines()
            .enumerate()
            .map(|(idx, line)| (idx.saturating_add(1), line.trim_end_matches('\r')))
            .filter(|(_, line)| !line.trim().is_empty());

        let (_, header) = lines
            .next()
            .ok_or_else(|| CuratomaticError::Red("the file is empty".to_owned()))?;
        let column = |name: &str| {
            header
                .split('\t')
                .position(|col| col.trim() == name)
                .ok_or_else(|| CuratomaticError::Red(format!("missing column: {name}")))
        };
        let rank_idx = column(COL_RANK)?;
        let median_idx = column(COL_MEDIAN)?;

        let mut data = BTreeMap::new();
        for (line_no, line) in lines {
            let cols: Vec<&str> = line.split('\t').collect();
            let value = |idx: usize, name: &str| {
                cols.get(idx).copied().ok_or_else(|| {
                    CuratomaticError::Red(format!("line {line_no} has no value for column: {name}"))
                })
            };
            let name = value(rank_idx, COL_RANK)?;
            let rank = Rank::from_name(name).map_err(|_unknown_rank| {
                CuratomaticError::Red(format!("line {line_no} has an unknown rank: {name}"))
            })?;
            let median = value(median_idx, COL_MEDIAN)?;
            let red: RedValue = median.parse().map_err(|e| {
                CuratomaticError::Red(format!(
                    "line {line_no} has an invalid median for {}: {median} ({e})",
                    rank.name()
                ))
            })?;

            if let Some(fixed) = fixed_red(rank)
                && red != fixed
            {
                return Err(CuratomaticError::Red(format!(
                    "line {line_no}: the {} median must be {fixed}, not {median}",
                    rank.name()
                )));
            }
            if data.insert(rank, red).is_some() {
                return Err(CuratomaticError::Red(format!(
                    "line {line_no}: {} is given more than once",
                    rank.name()
                )));
            }
        }

        for rank in [Rank::Domain, Rank::Species] {
            if let Some(fixed) = fixed_red(rank) {
                data.insert(rank, fixed);
            }
        }
        Ok(Red { data })
    }

    /// The RED value of a rank.
    pub(crate) fn get(&self, rank: Rank) -> CuratomaticResult<RedValue> {
        self.data
            .get(&rank)
            .copied()
            .ok_or(CuratomaticError::MissingRed(rank.prefix()))
    }

    /// Format as `<rank>=<red>`, ordered from the highest rank to the lowest.
    pub(crate) fn to_str(&self) -> String {
        self.data
            .iter()
            .map(|(rank, red)| format!("{}={red:.3}", rank.prefix()))
            .collect::<Vec<_>>()
            .join(", ")
    }
}

#[cfg(test)]
mod tests;
