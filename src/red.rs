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
        RedValue((self.0 + other.0) / 2.0)
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

/// The median relative evolutionary divergence (RED) expected for each rank.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Red {
    data: BTreeMap<Rank, RedValue>,
}

impl Red {
    /// Read a JSON object that maps a rank name to its RED value, e.g. `{"phylum": 0.23}`.
    ///
    /// The rank is determined by the first character of the name. The domain is always `0.0`.
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        let text = fs::read_to_string(path)
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
        Self::from_json(&text).map_err(|e| match e {
            CuratomaticError::Red(message) => {
                CuratomaticError::Red(format!("{}: {message}", path.display()))
            }
            other => other,
        })
    }

    /// Parse the JSON text that [`Red::from_file`] reads.
    pub(crate) fn from_json(text: &str) -> CuratomaticResult<Self> {
        let values: BTreeMap<String, f64> =
            serde_json::from_str(text).map_err(|e| CuratomaticError::Red(e.to_string()))?;

        let mut data = BTreeMap::new();
        for (name, value) in values {
            data.insert(Rank::of_taxon(&name)?, RedValue::new(value));
        }
        data.insert(Rank::Domain, RedValue::new(0.0));
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
