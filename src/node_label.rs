use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::newick::{NewickTree, NodeId};
use crate::red::RedValue;
use crate::taxonomy::Taxonomy;

/// Separates the rest of a node label from its RED value when it is written.
pub(crate) const RED_SEPARATOR: &str = "|RED=";

/// The spellings of the RED key that are read, which follow a `|` or a `:`.
const RED_KEYS: [&str; 2] = ["RED=", "red="];

/// The information held in the label of a node.
///
/// A leaf is labelled `<gid>[|RED=<red>]`, a leaf without a RED has a RED of [`RedValue::LEAF`].
/// A leaf label never holds a taxonomy, and is never written: the taxa of a leaf are kept on
/// [`crate::tree::Leaf`]. An internal node is labelled `[<bootstrap>][:][<taxonomy>]|RED=<red>`, with
/// at least one of the bootstrap and taxonomy. The taxonomy holds named taxa only: a placeholder
/// (e.g. `100:g__|RED=0.9`) is [`CuratomaticError::PlaceholderOnTree`]. The RED key may also be written `|red=`, or `:red=`
/// straight after a bootstrap (e.g. `100:red=0.889`).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NodeLabel {
    pub(crate) red: RedValue,
    pub(crate) bootstrap: Option<BootstrapValue>,
    pub(crate) taxonomy: Option<Taxonomy>,
    pub(crate) gid: Option<GenomeId>,
}

impl NodeLabel {
    /// Read the label of a node. The `taxonomy` of a leaf is always `None`, its taxa are kept on
    /// [`crate::tree::Leaf`].
    pub(crate) fn from_node(tree: &NewickTree, node: NodeId) -> CuratomaticResult<Self> {
        let label = tree
            .label(node)
            .ok_or_else(|| invalid(&format!("node {node} has no label")))?;
        if tree.is_leaf(node) {
            return parse_leaf_label(label);
        }

        // Either the bootstrap, the taxonomy, or both.
        let (red, bootstrap, taxonomy) = parse_node_label(label)?;
        let taxonomy = parse_taxonomy(taxonomy.unwrap_or_default()).map_err(|e| {
            if matches!(e, CuratomaticError::Placeholder(_)) {
                CuratomaticError::PlaceholderOnTree(label.to_owned())
            } else {
                e
            }
        })?;
        Ok(NodeLabel {
            red,
            bootstrap,
            taxonomy,
            gid: None,
        })
    }

    /// Format as the label of an internal node, the inverse of [`NodeLabel::from_node`].
    ///
    /// A leaf label is never written (it is always the genome ID it was read with), so formatting
    /// the label of a leaf is an error.
    pub(crate) fn to_label(&self) -> CuratomaticResult<String> {
        if let Some(gid) = &self.gid {
            return Err(CuratomaticError::Invariant(format!(
                "the label of leaf {gid} is never rewritten"
            )));
        }
        let tax_string = self.taxonomy.as_ref().map_or_default(ToString::to_string);
        let red = self.red;
        match (self.bootstrap, tax_string.is_empty()) {
            (None, true) => Err(invalid("node has neither a bootstrap nor a taxonomy")),
            (None, false) => Ok(format!("{tax_string}{RED_SEPARATOR}{red:.3}")),
            (Some(bootstrap), true) => Ok(format!("{bootstrap:.1}{RED_SEPARATOR}{red:.3}")),
            (Some(bootstrap), false) => Ok(format!(
                "{bootstrap:.1}:{tax_string}{RED_SEPARATOR}{red:.3}"
            )),
        }
    }
}

/// Read the label of a leaf: `<gid>[|RED=<red>]`. Anything else after the genome ID (e.g. a
/// taxonomy, `<gid>|g__A`) is an error.
fn parse_leaf_label(label: &str) -> CuratomaticResult<NodeLabel> {
    let (gid, red) = split_red(label)?;
    let red = red.unwrap_or(RedValue::LEAF);
    if red != RedValue::LEAF || gid.is_empty() {
        return Err(invalid(label));
    }
    if gid.contains('|') {
        return Err(CuratomaticError::TaxonomyOnLeaf(label.to_owned()));
    }
    Ok(NodeLabel {
        red,
        bootstrap: None,
        taxonomy: None,
        gid: Some(GenomeId::new(gid)),
    })
}

/// `None` for an empty taxonomy string.
fn parse_taxonomy(taxonomy: &str) -> CuratomaticResult<Option<Taxonomy>> {
    if taxonomy.trim().is_empty() {
        Ok(None)
    } else {
        Taxonomy::from_string(taxonomy).map(Some)
    }
}

/// Parse an internal node label into its RED, bootstrap and taxonomy (if they are set).
fn parse_node_label(
    label: &str,
) -> CuratomaticResult<(RedValue, Option<BootstrapValue>, Option<&str>)> {
    let float = |value: &str| {
        value
            .parse::<BootstrapValue>()
            .map_err(|_not_a_float| invalid(label))
    };

    let (prefix, red) = split_red(label)?;
    let red = red.ok_or_else(|| invalid(label))?;

    // Either the taxonomy, the bootstrap, or the bootstrap and taxonomy.
    if prefix.contains(':') {
        let mut parts = prefix.split(':');
        let (Some(bootstrap), Some(taxon), None) = (parts.next(), parts.next(), parts.next())
        else {
            return Err(invalid(label));
        };
        Ok((red, Some(float(bootstrap)?), Some(taxon)))
    } else if let Ok(bootstrap) = float(prefix) {
        Ok((red, Some(bootstrap), None))
    } else {
        Ok((red, None, Some(prefix)))
    }
}

/// Split a label into what precedes the RED, and the RED value if there is one.
///
/// The RED follows `|RED=` or `|red=`, and the `|` is not part of the prefix. It can also follow
/// `:RED=` or `:red=`, where the `:` separates a bootstrap from an empty taxonomy, and is kept in
/// the prefix (`100:red=0.5` gives `100:`). A label with more than one RED is invalid.
fn split_red(label: &str) -> CuratomaticResult<(&str, Option<RedValue>)> {
    let label = label.trim();
    let mut found = RED_KEYS
        .iter()
        .flat_map(|key| label.match_indices(key))
        .map(|(idx, key)| (idx, key.len()));
    let Some((idx, key_len)) = found.next() else {
        return Ok((label, None));
    };
    if found.next().is_some() {
        return Err(invalid(label));
    }

    let (before, rest) = label.split_at_checked(idx).ok_or_else(|| invalid(label))?;
    let value = rest.get(key_len..).ok_or_else(|| invalid(label))?;
    let red = value.parse().map_err(|_not_a_float| invalid(label))?;
    let prefix = if let Some(prefix) = before.strip_suffix('|') {
        prefix
    } else if before.ends_with(':') {
        before
    } else {
        return Err(invalid(label));
    };
    Ok((prefix, Some(red)))
}

/// The label without its RED value, e.g. `100:g__A|RED=0.9` gives `100:g__A`, and `100:red=0.9`
/// gives `100`.
pub(crate) fn strip_red(label: &str) -> CuratomaticResult<&str> {
    let (prefix, red) = split_red(label)?;
    Ok(if red.is_some() {
        prefix.strip_suffix(':').unwrap_or(prefix)
    } else {
        prefix
    })
}

fn invalid(label: &str) -> CuratomaticError {
    CuratomaticError::InvalidNodeLabel(label.to_owned())
}

#[cfg(test)]
mod tests;
