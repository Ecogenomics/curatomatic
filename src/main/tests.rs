use super::*;

// Python formats counts with `f'{n:,}'`.
#[test]
fn thousands_matches_python_comma_formatting() {
    let cases = [
        (0, "0"),
        (7, "7"),
        (999, "999"),
        (1000, "1,000"),
        (12_345, "12,345"),
        (100_000, "100,000"),
        (1_234_567, "1,234,567"),
    ];
    for (n, expected) in cases {
        assert_eq!(thousands(n), expected, "n = {n}");
    }
}
