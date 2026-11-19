use super::*;

/// Parse a tree that is expected to be valid.
fn parse(text: &str) -> NewickTree {
    NewickTree::parse(text).unwrap()
}

// The expected values below were produced by dendropy 5.0.8 (`escape_nexus_token` with the
// arguments used by `dendropy/dataio/newickwriter.py::_render_node_tag`).
#[test]
fn escape_label_matches_dendropy() {
    let cases = [
        ("D-G1", "D-G1"),
        ("G1", "G1"),
        ("100.0", "100.0"),
        ("a=b", "a=b"),
        ("x y", "x_y"),
        // A tab is protected, so unlike a space it is quoted rather than replaced.
        ("x\ty", "'x\ty'"),
        ("s__A b", "'s__A b'"),
        ("g__X", "'g__X'"),
        ("a,b", "'a,b'"),
        ("it's", "'it''s'"),
        ("", ""),
    ];
    for (label, expected) in cases {
        assert_eq!(escape_label(label), expected, "label {label:?}");
    }
}

#[test]
fn parse_reads_quoted_labels_with_doubled_quotes() {
    let tree = parse("('it''s':1.0,'a b':2.0)'x';");
    let leaves: Vec<_> = tree.leaves().collect();
    assert_eq!(tree.label(leaves[0]), Some("it's"));
    assert_eq!(tree.label(leaves[1]), Some("a b"));
    assert_eq!(tree.label(tree.root), Some("x"));
}

#[test]
fn parse_keeps_underscores_in_unquoted_labels() {
    let tree = parse("(a_b:1.0,c:2.0);");
    let leaves: Vec<_> = tree.leaves().collect();
    assert_eq!(tree.label(leaves[0]), Some("a_b"));
}

#[test]
fn parse_skips_comments() {
    let tree = parse("[&R] (A:1.0,[a comment]B:2.0)[&x=1]C;");
    let leaves: Vec<_> = tree.leaves().collect();
    assert_eq!(leaves.len(), 2);
    assert_eq!(tree.label(leaves[1]), Some("B"));
    assert_eq!(tree.label(tree.root), Some("C"));
}

#[test]
fn parse_rejects_unbalanced_parentheses() {
    assert!(matches!(
        NewickTree::parse("((A,B);"),
        Err(CuratomaticError::Newick(_))
    ));
    assert!(matches!(
        NewickTree::parse("(A,B));"),
        Err(CuratomaticError::Newick(_))
    ));
}

#[test]
fn parse_rejects_bad_branch_length() {
    assert!(matches!(
        NewickTree::parse("(A:x,B);"),
        Err(CuratomaticError::Newick(_))
    ));
}

#[test]
fn parse_rejects_unterminated_quote() {
    assert!(matches!(
        NewickTree::parse("('A,B);"),
        Err(CuratomaticError::Newick(_))
    ));
}

#[test]
fn ids_are_post_order_with_the_root_last() {
    // ((A,B)X,(C,D)Y)R, in post-order: A B X C D Y R.
    let tree = parse("((A,B)X,(C,D)Y)R;");
    let labels: Vec<_> = tree.ids().map(|id| tree.label(id).unwrap()).collect();
    assert_eq!(labels, ["A", "B", "X", "C", "D", "Y", "R"]);
    assert_eq!(tree.ids().last(), Some(tree.root));
}

#[test]
fn leaves_are_left_to_right() {
    let tree = parse("((A,B)X,(C,(D,E)Z)Y)R;");
    let labels: Vec<_> = tree.leaves().map(|id| tree.label(id).unwrap()).collect();
    assert_eq!(labels, ["A", "B", "C", "D", "E"]);
}

#[test]
fn descendants_are_the_contiguous_subtree_range() {
    let tree = parse("((A,B)X,(C,(D,E)Z)Y)R;");
    let find = |label: &str| {
        tree.ids()
            .find(|&id| tree.label(id) == Some(label))
            .unwrap()
    };
    let labels = |id| -> Vec<&str> {
        tree.descendants(id)
            .map(|d| tree.label(d).unwrap())
            .collect()
    };
    assert_eq!(labels(find("A")), Vec::<&str>::new());
    assert_eq!(labels(find("X")), ["A", "B"]);
    assert_eq!(labels(find("Y")), ["C", "D", "E", "Z"]);
    assert_eq!(labels(find("R")).len(), 8);
}

#[test]
fn parent_and_children_agree() {
    let tree = parse("((A,B)X,(C,(D,E)Z)Y)R;");
    assert_eq!(tree.parent(tree.root), None);
    for id in tree.ids() {
        for &child in tree.children(id) {
            assert_eq!(tree.parent(child), Some(id));
        }
        assert_eq!(tree.is_leaf(id), tree.children(id).is_empty());
    }
    let root_children: Vec<_> = tree
        .children(tree.root)
        .iter()
        .map(|&c| tree.label(c).unwrap())
        .collect();
    assert_eq!(root_children, ["X", "Y"]);
}

#[test]
fn set_label_replaces_the_label() {
    let mut tree = parse("(A,B)R;");
    let root = tree.root;
    tree.set_label(root, "s__A b".to_string());
    assert_eq!(tree.to_newick(), "(A,B)'s__A b';\n");
}

// The input is already in the form dendropy writes, so reading and writing it is the identity.
// dendropy 5.0.8 gives the same string for the equivalent un-normalised input
// `((G1:0.10,D-G1:1)'s__A b':0.5,('xy':1e-5,'a=b':0.0,'it''s':3)'100.0':2)'d__X';`.
#[test]
fn round_trips_a_tree_in_dendropy_form() {
    let text = "((G1:0.1,D-G1:1.0)'s__A b':0.5,(xy:1e-05,a=b:0.0,'it''s':3.0)100.0:2.0)'d__X';";
    assert_eq!(parse(text).to_newick(), format!("{text}\n"));
}

#[test]
fn writing_normalises_label_quoting() {
    assert_eq!(
        parse("('x y','a=b','100.0')'d__X';").to_newick(),
        "(x_y,a=b,100.0)'d__X';\n"
    );
}

// dendropy 5.0.8, read with `preserve_underscores=True`: `(x_y:1.0,b:2.0)r;` is written as
// `('x_y':1.0,b:2.0)r;`.
#[test]
fn bare_underscore_labels_are_quoted_when_written() {
    assert_eq!(
        parse("(x_y:1.0,b:2.0)r;").to_newick(),
        "('x_y':1.0,b:2.0)r;\n"
    );
}

#[test]
fn node_ids_order_numerically_in_post_order_and_display_their_position() {
    // Post-order: A B X C D Y R, so ids 0..=6; 10 > 9 numerically, not as text.
    let tree = parse("((A,B)X,(C,D)Y)R;");
    let ids: Vec<NodeId> = tree.ids().collect();
    assert!(ids.windows(2).all(|pair| pair[0] < pair[1]));
    assert_eq!(
        ids.iter().map(|id| id.index()).collect::<Vec<_>>(),
        [0, 1, 2, 3, 4, 5, 6]
    );
    assert_eq!(ids.last().unwrap().to_string(), "6");
    assert_eq!(tree.root.to_string(), "6");

    let big = format!(
        "({});",
        (0..11)
            .map(|i| format!("L{i}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    let big = parse(&big);
    let ids: Vec<NodeId> = big.ids().collect();
    assert!(ids[9] < ids[10]);
    assert_eq!(ids[10].to_string(), "10");
}
