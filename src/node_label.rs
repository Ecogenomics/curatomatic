use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::GenomeId;
use crate::newick::{NewickTree, NodeId};
use crate::red::RedValue;
use crate::taxonomy::Taxonomy;
use crate::tree::parse_node_label;

/// Separates the rest of a node label from its RED value.
const RED_SEPARATOR: &str = "|RED=";

/// The information held in the label of a node.
///
/// A leaf is labelled `<gid>|RED=<red>`. An internal node is labelled
/// `[<bootstrap>][:][<taxonomy>]|RED=<red>`, with at least one of the bootstrap and taxonomy.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NodeLabel {
    pub(crate) red: RedValue,
    pub(crate) bootstrap: Option<BootstrapValue>,
    pub(crate) taxonomy: Option<Taxonomy>,
    pub(crate) gid: Option<GenomeId>,
}

impl NodeLabel {
    /// Read the label of a node.
    pub(crate) fn from_node(tree: &NewickTree, node: NodeId) -> CuratomaticResult<Self> {
        let label = tree
            .label(node)
            .ok_or_else(|| invalid(&format!("node {node} has no label")))?;
        if tree.is_leaf(node) {
            let (prefix, red) = split_red(label)?;
            if red != RedValue::LEAF {
                return Err(invalid(label));
            }
            return Ok(NodeLabel {
                red,
                bootstrap: None,
                taxonomy: None,
                gid: Some(GenomeId::new(prefix)),
            });
        }

        // Either the bootstrap, the taxonomy, or both.
        let (red, bootstrap, taxonomy) = parse_node_label(label)?;
        let taxonomy = taxonomy.unwrap_or_default();

        let taxonomy = if taxonomy.is_empty() {
            None
        } else {
            Some(Taxonomy::from_string(taxonomy)?)
        };
        Ok(NodeLabel {
            red,
            bootstrap,
            taxonomy,
            gid: None,
        })
    }

    /// Format as the label of a node, the inverse of [`NodeLabel::from_node`].
    pub(crate) fn to_label(&self) -> CuratomaticResult<String> {
        let tax_string = self.taxonomy.as_ref().map_or_default(ToString::to_string);
        let red = self.red;

        if let Some(gid) = &self.gid {
            // This is a leaf node.
            if !tax_string.is_empty() {
                return Err(invalid(&format!(
                    "leaf {gid} cannot hold a taxonomy: {tax_string}"
                )));
            }
            return Ok(format!("{gid}{RED_SEPARATOR}{red:.3}"));
        }

        // This is an internal node.
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

/// Split a label into what precedes `|RED=`, and the RED value.
fn split_red(label: &str) -> CuratomaticResult<(&str, RedValue)> {
    let mut parts = label.trim().split(RED_SEPARATOR);
    match (parts.next(), parts.next(), parts.next()) {
        (Some(prefix), Some(red), None) => Ok((prefix, parse_red(red, label)?)),
        _ => Err(invalid(label)),
    }
}

fn parse_red(value: &str, label: &str) -> CuratomaticResult<RedValue> {
    value.parse().map_err(|_not_a_float| invalid(label))
}

fn invalid(label: &str) -> CuratomaticError {
    CuratomaticError::InvalidNodeLabel(label.to_owned())
}

#[cfg(test)]
mod tests;
