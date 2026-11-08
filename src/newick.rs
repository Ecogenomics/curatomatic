//! A minimal Newick tree, standing in for `dendropy.Tree` (read with `preserve_underscores=True`,
//! written with `suppress_rooting=True`).
//!
//! Nodes live in an arena and are numbered in post-order: every node comes after all of its
//! descendants, and the descendants of a node are the contiguous ids directly before it. Iterating
//! over ids therefore matches `postorder_node_iter`, and iterating over the leaves matches
//! `leaf_node_iter`.

use std::borrow::Cow;
use std::fmt;
use std::fs;
use std::path::Path;

use crate::error::{CuratomaticError, CuratomaticResult};

/// Identifies a node of a [`NewickTree`].
///
/// Ids are numbered in post-order, so the derived order is the post-order of the tree. Only the
/// tree creates them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct NodeId(usize);

impl NodeId {
    /// The position of the node in post-order, to index data kept for every node.
    pub(crate) fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// Characters that force a label to be quoted when written: the `protect_regex` that dendropy's
/// Newick writer passes to `escape_nexus_token`, `[()[\]{},;:'"\0\t\n]`.
const PROTECTED_CHARS: &str = "()[]{},;:'\"\0\t\n";

#[derive(Debug, Clone, PartialEq, Eq)]
struct Node {
    /// The label of a leaf (its taxon) or of an internal node.
    label: Option<String>,
    /// The branch length, exactly as it was read.
    length: Option<String>,
    parent: Option<NodeId>,
    children: Vec<NodeId>,
    /// The id of the first node in this node's subtree.
    first: NodeId,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct NewickTree {
    nodes: Vec<Node>,
    root: NodeId,
}

impl NewickTree {
    /// Read a Newick file.
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        let text = fs::read_to_string(path)
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
        Self::parse(&text)
    }

    /// Parse a Newick string. Comments (`[...]`, e.g. `[&R]`) are ignored.
    pub(crate) fn parse(text: &str) -> CuratomaticResult<Self> {
        let mut parser = Parser {
            chars: text.chars().collect(),
            pos: 0,
        };
        let mut nodes: Vec<Node> = Vec::new();
        // The children collected so far for each parenthesis that is still open.
        let mut open: Vec<Vec<NodeId>> = Vec::new();

        let root = 'outer: loop {
            parser.skip_blank();
            if parser.eat('(') {
                open.push(Vec::new());
                continue;
            }

            let (label, length) = parser.read_label_and_length()?;
            let mut id = NodeId(nodes.len());
            nodes.push(Node {
                label,
                length,
                parent: None,
                children: Vec::new(),
                first: id,
            });

            // Attach the finished node, closing every parenthesis that ends after it.
            loop {
                let Some(siblings) = open.last_mut() else {
                    break 'outer id;
                };
                siblings.push(id);

                parser.skip_blank();
                if parser.eat(',') {
                    continue 'outer;
                }
                if !parser.eat(')') {
                    return Err(parser.error("expected `,` or `)`"));
                }

                let children = open.pop().unwrap_or_default();
                let (label, length) = parser.read_label_and_length()?;
                id = NodeId(nodes.len());
                let first = children
                    .first()
                    .and_then(|&child| nodes.get(child.0))
                    .map_or(id, |child| child.first);
                for &child in &children {
                    if let Some(node) = nodes.get_mut(child.0) {
                        node.parent = Some(id);
                    }
                }
                nodes.push(Node {
                    label,
                    length,
                    parent: None,
                    children,
                    first,
                });
            }
        };

        parser.skip_blank();
        parser.eat(';');
        parser.skip_blank();
        if parser.pos < parser.chars.len() {
            return Err(parser.error("unexpected trailing characters"));
        }

        Ok(NewickTree { nodes, root })
    }

    /// Format as a Newick string, terminated by `;` and a newline.
    pub(crate) fn to_newick(&self) -> String {
        enum Step {
            Enter(NodeId),
            Separator,
            Close(NodeId),
        }

        let mut out = String::new();
        let mut stack = vec![Step::Enter(self.root)];
        while let Some(step) = stack.pop() {
            match step {
                Step::Enter(id) => {
                    let Some(node) = self.nodes.get(id.0) else {
                        continue;
                    };
                    if node.children.is_empty() {
                        write_label_and_length(&mut out, node);
                    } else {
                        out.push('(');
                        stack.push(Step::Close(id));
                        for (i, &child) in node.children.iter().enumerate().rev() {
                            stack.push(Step::Enter(child));
                            if i > 0 {
                                stack.push(Step::Separator);
                            }
                        }
                    }
                }
                Step::Separator => out.push(','),
                Step::Close(id) => {
                    out.push(')');
                    if let Some(node) = self.nodes.get(id.0) {
                        write_label_and_length(&mut out, node);
                    }
                }
            }
        }
        out.push_str(";\n");
        out
    }

    /// Write the tree to a file, see [`NewickTree::to_newick`].
    pub(crate) fn to_file(&self, path: &Path) -> CuratomaticResult<()> {
        fs::write(path, self.to_newick())
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))
    }

    /// Every node, in post-order.
    pub(crate) fn ids(&self) -> impl Iterator<Item = NodeId> + use<> {
        (0..self.nodes.len()).map(NodeId)
    }

    /// Every leaf, from left to right.
    pub(crate) fn leaves(&self) -> impl Iterator<Item = NodeId> {
        self.ids().filter(|&id| self.is_leaf(id))
    }

    pub(crate) fn is_leaf(&self, id: NodeId) -> bool {
        self.nodes.get(id.0).is_some_and(|node| node.children.is_empty())
    }

    pub(crate) fn parent(&self, id: NodeId) -> Option<NodeId> {
        self.nodes.get(id.0)?.parent
    }

    pub(crate) fn children(&self, id: NodeId) -> &[NodeId] {
        self.nodes
            .get(id.0)
            .map_or(&[], |node| node.children.as_slice())
    }

    /// Every node below this one, excluding the node itself.
    pub(crate) fn descendants(&self, id: NodeId) -> impl Iterator<Item = NodeId> + use<> {
        let range = self
            .nodes
            .get(id.0)
            .map_or(0..0, |node| node.first.0..id.0);
        range.map(NodeId)
    }

    pub(crate) fn label(&self, id: NodeId) -> Option<&str> {
        self.nodes.get(id.0)?.label.as_deref()
    }

    pub(crate) fn set_label(&mut self, id: NodeId, label: String) {
        if let Some(node) = self.nodes.get_mut(id.0) {
            node.label = Some(label);
        }
    }
}

/// Escape a label as dendropy's Newick writer does: `escape_nexus_token(label,
/// preserve_spaces=False, quote_underscores=True, protect_regex=...)` in
/// `dendropy/dataio/nexusprocessing.py`, with `protect_regex` as in [`PROTECTED_CHARS`].
///
/// A label with no `_` and no protected character is written bare, with spaces replaced by `_`.
/// Any other label is quoted, with `'` doubled.
fn escape_label(label: &str) -> Cow<'_, str> {
    if !label.contains('_') && !label.contains(|c| PROTECTED_CHARS.contains(c)) {
        if label.contains(' ') {
            Cow::Owned(label.replace(' ', "_"))
        } else {
            Cow::Borrowed(label)
        }
    } else {
        Cow::Owned(format!("'{}'", label.replace('\'', "''")))
    }
}

fn write_label_and_length(out: &mut String, node: &Node) {
    if let Some(label) = &node.label {
        out.push_str(&escape_label(label));
    }
    if let Some(length) = &node.length {
        out.push(':');
        out.push_str(length);
    }
}

struct Parser {
    chars: Vec<char>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn eat(&mut self, expected: char) -> bool {
        let found = self.peek() == Some(expected);
        if found {
            self.pos = self.pos.saturating_add(1);
        }
        found
    }

    fn error(&self, message: &str) -> CuratomaticError {
        CuratomaticError::Newick(format!("{message} (at character {})", self.pos))
    }

    /// Skip whitespace and comments.
    fn skip_blank(&mut self) {
        while let Some(c) = self.peek() {
            if c.is_whitespace() {
                self.pos = self.pos.saturating_add(1);
            } else if c == '[' {
                while let Some(c) = self.peek() {
                    self.pos = self.pos.saturating_add(1);
                    if c == ']' {
                        break;
                    }
                }
            } else {
                break;
            }
        }
    }

    /// Read an optional label followed by an optional `:<length>`.
    fn read_label_and_length(&mut self) -> CuratomaticResult<(Option<String>, Option<String>)> {
        self.skip_blank();
        let label = if self.eat('\'') {
            Some(self.read_quoted()?)
        } else {
            let label = self.read_while(|c| !c.is_whitespace() && !"(),:;[]'".contains(c));
            (!label.is_empty()).then_some(label)
        };

        self.skip_blank();
        let length = if self.eat(':') {
            self.skip_blank();
            let length = self.read_while(|c| !c.is_whitespace() && !"(),;[]".contains(c));
            if length.parse::<f64>().is_err() {
                return Err(self.error(&format!("invalid branch length `{length}`")));
            }
            Some(length)
        } else {
            None
        };
        Ok((label, length))
    }

    /// Read up to the closing quote, the opening one has been consumed. `''` is a literal quote.
    fn read_quoted(&mut self) -> CuratomaticResult<String> {
        let mut out = String::new();
        while let Some(c) = self.peek() {
            self.pos = self.pos.saturating_add(1);
            if c == '\'' {
                if self.eat('\'') {
                    out.push('\'');
                } else {
                    return Ok(out);
                }
            } else {
                out.push(c);
            }
        }
        Err(self.error("unterminated quoted label"))
    }

    fn read_while(&mut self, keep: impl Fn(char) -> bool) -> String {
        let mut out = String::new();
        while let Some(c) = self.peek()
            && keep(c)
        {
            out.push(c);
            self.pos = self.pos.saturating_add(1);
        }
        out
    }
}

#[cfg(test)]
mod tests;
