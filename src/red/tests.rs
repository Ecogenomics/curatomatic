use super::*;

const JSON: &str = r#"{"phylum": 0.3, "class": 0.5, "order": 0.7, "family": 0.8, "genus": 0.9}"#;

fn red() -> Red {
    Red::from_json(JSON).unwrap()
}

#[test]
fn ranks_are_keyed_on_the_first_letter() {
    // `curatomatic/model/red.py` uses `k[0]`, so `p` and `phylum` are the same key.
    let red = Red::from_json(r#"{"p": 0.1, "cxyz": 0.2}"#).unwrap();
    assert_eq!(red.get(Rank::Phylum), Ok(RedValue::new(0.1)));
    assert_eq!(red.get(Rank::Class), Ok(RedValue::new(0.2)));
}

#[test]
fn a_name_that_is_not_a_rank_is_an_error() {
    // Python stores this under the key `x`, where nothing ever reads it. Rust rejects it.
    assert_eq!(
        Red::from_json(r#"{"xyz": 0.2}"#),
        Err(CuratomaticError::UnknownRank("xyz".to_owned()))
    );
}

#[test]
fn the_domain_is_always_zero() {
    let red = red();
    assert_eq!(red.get(Rank::Domain), Ok(RedValue::new(0.0)));

    // A domain in the input is overwritten, like `out['d'] = 0.0` in Python (`model/red.py`).
    let red = Red::from_json(r#"{"domain": 0.4}"#).unwrap();
    assert_eq!(red.get(Rank::Domain), Ok(RedValue::new(0.0)));
}

#[test]
fn a_missing_rank_is_an_error() {
    let red = red();
    assert_eq!(red.get(Rank::Species), Err(CuratomaticError::MissingRed('s')));
}

// Expected value from Python: `RED({...}).to_str()` is
// `d=0.000, p=0.300, c=0.500, o=0.700, f=0.800, g=0.900`.
#[test]
fn to_str_is_in_rank_order_with_three_decimals() {
    assert_eq!(
        red().to_str(),
        "d=0.000, p=0.300, c=0.500, o=0.700, f=0.800, g=0.900"
    );

    // The order of the keys in the input does not matter.
    let red = Red::from_json(r#"{"genus": 0.9, "phylum": 0.3}"#).unwrap();
    assert_eq!(red.to_str(), "d=0.000, p=0.300, g=0.900");
}

#[test]
fn to_str_rounds_like_python() {
    // Python `'%.3f' % x`: 0.0625 -> 0.062, 0.0635 -> 0.064, 0.9995 -> 1.000.
    let red = Red::from_json(r#"{"phylum": 0.0625, "class": 0.0635, "order": 0.9995}"#).unwrap();
    assert_eq!(red.to_str(), "d=0.000, p=0.062, c=0.064, o=1.000");
}

#[test]
fn invalid_json_is_an_error() {
    assert!(matches!(Red::from_json("{"), Err(CuratomaticError::Red(_))));
    assert!(matches!(
        Red::from_json(r#"{"phylum": "x"}"#),
        Err(CuratomaticError::Red(_))
    ));
}

#[test]
fn an_empty_key_is_an_error() {
    assert_eq!(Red::from_json(r#"{"": 0.1}"#), Err(CuratomaticError::EmptyTaxon));
}

#[test]
fn from_file_reads_the_json() {
    let path = std::env::temp_dir().join(format!("curatomatic_red_{}.json", std::process::id()));
    fs::write(&path, JSON).unwrap();
    let from_file = Red::from_file(&path);
    fs::remove_file(&path).unwrap();
    assert_eq!(from_file.unwrap(), red());
}

#[test]
fn from_file_reports_the_path_of_a_bad_file() {
    let path = std::env::temp_dir().join(format!("curatomatic_red_bad_{}.json", std::process::id()));
    fs::write(&path, "{").unwrap();
    let from_file = Red::from_file(&path);
    fs::remove_file(&path).unwrap();
    match from_file {
        Err(CuratomaticError::Red(message)) => {
            assert!(message.contains(&path.display().to_string()), "{message}");
        }
        other => panic!("expected a Red error, got {other:?}"),
    }
}

#[test]
fn from_file_with_a_missing_file_is_an_io_error() {
    assert!(matches!(
        Red::from_file(Path::new("does/not/exist.json")),
        Err(CuratomaticError::Io(_))
    ));
}

#[test]
fn midpoint_is_halfway() {
    let low = RedValue::new(0.5);
    let high = RedValue::new(0.75);
    assert_eq!(low.midpoint(high), RedValue::new(0.625));
    assert_eq!(high.midpoint(low), RedValue::new(0.625));
    assert_eq!(low.midpoint(low), low);
}

#[test]
fn distance_is_absolute() {
    let low = RedValue::new(0.25);
    let high = RedValue::new(0.75);
    assert!((low.distance(high) - 0.5).abs() < f64::EPSILON);
    assert!((high.distance(low) - 0.5).abs() < f64::EPSILON);
    assert!(low.distance(low).abs() < f64::EPSILON);
}

#[test]
fn red_value_parses_and_trims() {
    assert_eq!(" 0.25 ".parse(), Ok(RedValue::new(0.25)));
    assert!("x".parse::<RedValue>().is_err());
    assert_eq!(RedValue::LEAF, RedValue::new(1.0));
}

#[test]
fn red_value_display_respects_the_precision() {
    assert_eq!(format!("{:.3}", RedValue::new(0.957_017)), "0.957");
}

#[test]
fn halfway_to_matches_the_python_arithmetic() {
    // Python: 0.098 + ((0.472 - 0.098) / 2) == 0.28500000000000003, but (0.098 + 0.472) / 2 == 0.285.
    let (parent, rank) = (RedValue::new(0.098), RedValue::new(0.472));
    assert_eq!(parent.halfway_to(rank), RedValue::new(0.285_000_000_000_000_03));
    assert_eq!(parent.midpoint(rank), RedValue::new(0.285));
    assert!(RedValue::new(0.285) < parent.halfway_to(rank));
    assert_eq!(parent.halfway_to(parent), parent);
}
