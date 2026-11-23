use std::cmp::Ordering;
use std::fmt;
use std::num::ParseFloatError;
use std::str::FromStr;

/// A bootstrap support value.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct BootstrapValue(f64);

impl BootstrapValue {
    /// The bootstrap support assumed for a leaf node.
    pub(crate) const LEAF: BootstrapValue = BootstrapValue(100.0);

    pub(crate) const fn new(value: f64) -> Self {
        BootstrapValue(value)
    }

    /// A total ordering of bootstrap values, see [`f64::total_cmp`].
    pub(crate) fn total_cmp(self, other: BootstrapValue) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

impl FromStr for BootstrapValue {
    type Err = ParseFloatError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        value.trim().parse().map(BootstrapValue)
    }
}

/// Formats like the underlying `f64`, so a precision such as `{:.1}` is respected.
impl fmt::Display for BootstrapValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, f)
    }
}

#[cfg(test)]
mod tests;
