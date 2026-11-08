use super::*;

// Expected values in this file come from the Python implementation
// (`curatomatic/model/node_label.py`, `NodeLabel.from_node` and `to_string`), run with dendropy
// on the same labels.

/// The label of the root of a two leaf tree whose root is labelled `label`.
fn internal(label: &str) -> CuratomaticResult<NodeLabel> {
    let tree = NewickTree::parse(&format!(
        "('G1|RED=1.0':1.0,'G2|RED=1.0':1.0)'{label}';"
    ))
    .unwrap();
    // The root is the last node in post-order.
    let root = tree.ids().last().unwrap();
    NodeLabel::from_node(&tree, root)
}

/// The label of the first leaf of a two leaf tree whose first leaf is labelled `label`.
fn leaf(label: &str) -> CuratomaticResult<NodeLabel> {
    let tree = NewickTree::parse(&format!(
        "('{label}':1.0,'G2|RED=1.0':1.0)'50.0|RED=0.5';"
    ))
    .unwrap();
    let first = tree.leaves().next().unwrap();
    NodeLabel::from_node(&tree, first)
}

/// Read an internal node label, and write it back out.
fn round_trip(label: &str) -> String {
    internal(label).unwrap().to_label().unwrap()
}

#[test]
fn leaf_has_a_gid_and_a_red_of_one() {
    let label = leaf("G1|RED=1.0").unwrap();
    assert_eq!(
        label,
        NodeLabel {
            red: RedValue::LEAF,
            bootstrap: None,
            taxonomy: None,
            gid: Some(GenomeId::new("G1")),
        }
    );
    assert_eq!(label.to_label().unwrap(), "G1|RED=1.000");
}

#[test]
fn leaf_red_may_be_written_as_an_integer() {
    assert_eq!(leaf("G1|RED=1").unwrap().to_label().unwrap(), "G1|RED=1.000");
}

#[test]
fn leaf_with_a_red_other_than_one_is_an_error() {
    // Python: `AssertionError`.
    assert_eq!(
        leaf("G1|RED=0.5"),
        Err(CuratomaticError::InvalidNodeLabel("G1|RED=0.5".to_owned()))
    );
}

#[test]
fn leaf_without_a_red_or_with_a_bad_red_is_an_error() {
    // Python: `AssertionError` for the first, `ValueError` for the others.
    for label in ["G1", "G1|RED=x", "G1|RED=1.0|RED=1.0"] {
        assert_eq!(
            leaf(label),
            Err(CuratomaticError::InvalidNodeLabel(label.to_owned())),
            "label {label:?}"
        );
    }
}

#[test]
fn bootstrap_only() {
    let label = internal("100.0|RED=0.5").unwrap();
    assert_eq!(label.red, RedValue::new(0.5));
    assert_eq!(label.bootstrap, Some(BootstrapValue::new(100.0)));
    assert_eq!(label.taxonomy, None);
    assert_eq!(label.gid, None);
    assert_eq!(label.to_label().unwrap(), "100.0|RED=0.500");
}

#[test]
fn bootstrap_may_be_written_as_an_integer() {
    assert_eq!(round_trip("90|RED=0.25"), "90.0|RED=0.250");
}

#[test]
fn taxonomy_only() {
    let label = internal("f__X; g__Y|RED=0.5").unwrap();
    assert_eq!(label.red, RedValue::new(0.5));
    assert_eq!(label.bootstrap, None);
    assert_eq!(label.taxonomy, Some(Taxonomy::from_string("f__X; g__Y").unwrap()));
    assert_eq!(label.gid, None);
    assert_eq!(label.to_label().unwrap(), "f__X; g__Y|RED=0.500");
}

#[test]
fn bootstrap_and_taxonomy() {
    let label = internal("99.0:f__X|RED=0.5").unwrap();
    assert_eq!(label.red, RedValue::new(0.5));
    assert_eq!(label.bootstrap, Some(BootstrapValue::new(99.0)));
    assert_eq!(label.taxonomy, Some(Taxonomy::from_string("f__X").unwrap()));
    assert_eq!(label.gid, None);
    assert_eq!(label.to_label().unwrap(), "99.0:f__X|RED=0.500");
}

#[test]
fn a_red_of_one_on_an_internal_node_is_allowed() {
    assert_eq!(round_trip("100.0|RED=1.0"), "100.0|RED=1.000");
    assert_eq!(round_trip("f__X|RED=1"), "f__X|RED=1.000");
}

#[test]
fn surrounding_whitespace_is_ignored() {
    assert_eq!(round_trip(" 95.0|RED=0.5"), "95.0|RED=0.500");
}

#[test]
fn red_is_formatted_to_three_decimals_like_python() {
    // Python `'%.3f' % x`, an exact tie rounds to even.
    assert_eq!(round_trip("95.0|RED=0.0625"), "95.0|RED=0.062");
    assert_eq!(round_trip("95.0|RED=0.0635"), "95.0|RED=0.064");
    assert_eq!(round_trip("95.0|RED=0.957017"), "95.0|RED=0.957");
    assert_eq!(round_trip("95.0|RED=0.9995"), "95.0|RED=1.000");
    assert_eq!(round_trip("95.0|RED=0.0005"), "95.0|RED=0.001");
    assert_eq!(round_trip("95.0|RED=0.0015"), "95.0|RED=0.002");
    assert_eq!(round_trip("95.0|RED=0.0025"), "95.0|RED=0.003");
}

#[test]
fn bootstrap_is_formatted_to_one_decimal_like_python() {
    // Python `f'{x:.1f}'`.
    assert_eq!(round_trip("99.95|RED=0.5"), "100.0|RED=0.500");
    assert_eq!(round_trip("87.25|RED=0.5"), "87.2|RED=0.500");
}

#[test]
fn malformed_internal_labels_are_errors() {
    // Python: `AssertionError` for `noRED`, otherwise `ValueError`.
    for label in [
        "noRED",
        "x|RED=y",
        "a:b:c|RED=0.5",
        "abc:f__X|RED=0.5",
        "50|RED=0.5|RED=0.5",
    ] {
        assert_eq!(
            internal(label),
            Err(CuratomaticError::InvalidNodeLabel(label.to_owned())),
            "label {label:?}"
        );
    }
}

#[test]
fn internal_node_without_a_label_is_an_error() {
    let tree = NewickTree::parse("(A,B);").unwrap();
    let root = tree.ids().last().unwrap();
    assert!(matches!(
        NodeLabel::from_node(&tree, root),
        Err(CuratomaticError::InvalidNodeLabel(_))
    ));
}

#[test]
fn a_taxon_with_an_unknown_rank_is_an_error() {
    assert_eq!(
        internal("x__Foo|RED=0.5"),
        Err(CuratomaticError::UnknownRank("x__Foo".to_owned()))
    );
}

#[test]
fn to_label_of_a_node_with_neither_bootstrap_nor_taxonomy_is_an_error() {
    // Python: `AssertionError`.
    let label = NodeLabel {
        red: RedValue::new(0.5),
        bootstrap: None,
        taxonomy: None,
        gid: None,
    };
    assert!(matches!(
        label.to_label(),
        Err(CuratomaticError::InvalidNodeLabel(_))
    ));
}

#[test]
fn to_label_of_a_leaf_with_a_taxonomy_is_an_error() {
    // Python: `AssertionError`.
    let label = NodeLabel {
        red: RedValue::LEAF,
        bootstrap: None,
        taxonomy: Some(Taxonomy::from_string("g__Foo").unwrap()),
        gid: Some(GenomeId::new("G1")),
    };
    assert!(matches!(
        label.to_label(),
        Err(CuratomaticError::InvalidNodeLabel(_))
    ));
}
