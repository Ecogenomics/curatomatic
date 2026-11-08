use super::*;

#[test]
fn parses_and_trims_whitespace() {
    assert_eq!("95.0".parse(), Ok(BootstrapValue::new(95.0)));
    assert_eq!("  87.5\t".parse(), Ok(BootstrapValue::new(87.5)));
    assert_eq!("100".parse(), Ok(BootstrapValue::new(100.0)));
}

#[test]
fn rejects_text_that_is_not_a_number() {
    assert!("".parse::<BootstrapValue>().is_err());
    assert!("   ".parse::<BootstrapValue>().is_err());
    assert!("f__Foo".parse::<BootstrapValue>().is_err());
    assert!("95.0:x".parse::<BootstrapValue>().is_err());
}

// Python accepts the same non-finite spellings with `float()`.
#[test]
fn accepts_the_non_finite_spellings_of_f64() {
    assert!("inf".parse::<BootstrapValue>().is_ok());
    assert!("nan".parse::<BootstrapValue>().is_ok());
}

#[test]
fn display_respects_the_precision() {
    // Python `f'{x:.1f}'`: 99.95 -> '100.0', 87.25 -> '87.2' (ties are to even), 45.0 -> '45.0'.
    assert_eq!(format!("{:.1}", BootstrapValue::new(99.95)), "100.0");
    assert_eq!(format!("{:.1}", BootstrapValue::new(87.25)), "87.2");
    assert_eq!(format!("{:.1}", BootstrapValue::new(87.35)), "87.3");
    assert_eq!(format!("{:.1}", BootstrapValue::new(45.0)), "45.0");
    assert_eq!(format!("{}", BootstrapValue::new(95.0)), "95");
}

#[test]
fn total_cmp_orders_by_value() {
    let low = BootstrapValue::new(45.0);
    let high = BootstrapValue::new(95.0);
    assert_eq!(low.total_cmp(high), Ordering::Less);
    assert_eq!(high.total_cmp(low), Ordering::Greater);
    assert_eq!(high.total_cmp(high), Ordering::Equal);
}

#[test]
fn total_cmp_is_total_for_nan_and_signed_zero() {
    let nan = BootstrapValue::new(f64::NAN);
    assert_eq!(nan.total_cmp(nan), Ordering::Equal);
    assert_eq!(nan.total_cmp(BootstrapValue::new(f64::INFINITY)), Ordering::Greater);
    assert_eq!(
        BootstrapValue::new(-0.0).total_cmp(BootstrapValue::new(0.0)),
        Ordering::Less
    );
}

#[test]
fn partial_ord_compares_by_value() {
    assert!(BootstrapValue::new(80.0) < BootstrapValue::new(95.0));
    assert!(BootstrapValue::new(95.0) >= BootstrapValue::new(95.0));
}

#[test]
fn leaf_bootstrap_is_one_hundred() {
    assert_eq!(BootstrapValue::LEAF, BootstrapValue::new(100.0));
}
