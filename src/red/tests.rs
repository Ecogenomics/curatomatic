use strum::IntoEnumIterator;

use super::*;

const TSV: &str = "rank\tmedian\n\
                   Phylum\t0.3\n\
                   Class\t0.5\n\
                   Order\t0.7\n\
                   Family\t0.8\n\
                   Genus\t0.9\n";

fn red() -> Red {
    Red::from_tsv(TSV).unwrap()
}

/// A RED error, or a panic naming what was returned instead.
fn red_error(result: CuratomaticResult<Red>) -> String {
    match result {
        Err(CuratomaticError::Red(message)) => message,
        other => panic!("expected a Red error, got {other:?}"),
    }
}

#[test]
fn reads_the_median_column_by_name_and_ignores_the_others() {
    // The header and rows of `test_data/r237/median_red_rank_distribution.tsv`, trimmed.
    let text = "rank\tnum_taxa\tnum_taxa_pass_filter\tmedian\tmean\n\
                Phylum\t161\t46\t0.3095\t0.3314\n\
                Genus\t20066\t19408\t0.9189\t0.9175\n";
    let red = Red::from_tsv(text).unwrap();
    assert_eq!(red.get(Rank::Phylum), Ok(RedValue::new(0.3095)));
    assert_eq!(red.get(Rank::Genus), Ok(RedValue::new(0.9189)));

    // The columns can be in any order.
    let red = Red::from_tsv("median\tx\trank\n0.25\tfoo\tClass\n").unwrap();
    assert_eq!(red.get(Rank::Class), Ok(RedValue::new(0.25)));
}

#[test]
fn rank_names_ignore_case() {
    let red = Red::from_tsv("rank\tmedian\nphylum\t0.1\nCLASS\t0.2\nOrDeR\t0.3\n").unwrap();
    assert_eq!(red.get(Rank::Phylum), Ok(RedValue::new(0.1)));
    assert_eq!(red.get(Rank::Class), Ok(RedValue::new(0.2)));
    assert_eq!(red.get(Rank::Order), Ok(RedValue::new(0.3)));
}

#[test]
fn a_name_that_is_not_a_rank_is_an_error() {
    // The old JSON used only the first character, so `p` or `cxyz` were ranks. They are not now.
    for name in ["xyz", "p", "Phyla", ""] {
        let message = red_error(Red::from_tsv(&format!("rank\tmedian\n{name}\t0.2\n")));
        assert!(message.contains("unknown rank"), "{message}");
        assert!(message.contains("line 2"), "{message}");
    }
}

#[test]
fn the_domain_is_zero_and_the_species_is_a_leaf() {
    let red = red();
    assert_eq!(red.get(Rank::Domain), Ok(RedValue::new(0.0)));
    assert_eq!(red.get(Rank::Species), Ok(RedValue::LEAF));

    // Rows with those values are accepted.
    let red = Red::from_tsv("rank\tmedian\nDomain\t0\nSpecies\t1.0\n").unwrap();
    assert_eq!(red.get(Rank::Domain), Ok(RedValue::new(0.0)));
    assert_eq!(red.get(Rank::Species), Ok(RedValue::LEAF));
}

#[test]
fn a_domain_or_species_with_another_value_is_an_error() {
    for row in ["Domain\t0.4", "Species\t0.98"] {
        let message = red_error(Red::from_tsv(&format!("rank\tmedian\n{row}\n")));
        assert!(message.contains("must be"), "{message}");
    }
}

#[test]
fn a_rank_without_a_row_has_no_value() {
    let red = Red::from_tsv("rank\tmedian\nGenus\t0.9\n").unwrap();
    assert_eq!(red.get(Rank::Family), Err(CuratomaticError::MissingRed('f')));
}

#[test]
fn a_missing_column_is_an_error() {
    let message = red_error(Red::from_tsv("rank\tmean\nPhylum\t0.3\n"));
    assert!(message.contains("missing column: median"), "{message}");

    let message = red_error(Red::from_tsv("taxon\tmedian\nPhylum\t0.3\n"));
    assert!(message.contains("missing column: rank"), "{message}");

    let message = red_error(Red::from_tsv(""));
    assert!(message.contains("empty"), "{message}");
}

#[test]
fn a_row_too_short_for_a_column_is_an_error() {
    let message = red_error(Red::from_tsv("rank\tnum_taxa\tmedian\nPhylum\t161\n"));
    assert!(message.contains("line 2 has no value for column: median"), "{message}");
}

#[test]
fn a_median_that_is_not_a_number_is_an_error() {
    for median in ["x", "", "0.3.1"] {
        let message = red_error(Red::from_tsv(&format!("rank\tmedian\nPhylum\t{median}\n")));
        assert!(message.contains("line 2 has an invalid median for Phylum"), "{message}");
    }
}

#[test]
fn a_duplicate_rank_is_an_error() {
    // Also when the two rows differ only in case.
    let message = red_error(Red::from_tsv("rank\tmedian\nPhylum\t0.3\nGenus\t0.9\nphylum\t0.3\n"));
    assert!(message.contains("line 4: Phylum is given more than once"), "{message}");
}

#[test]
fn blank_lines_and_crlf_are_ignored() {
    let red = Red::from_tsv("rank\tmedian\r\n\r\nPhylum\t0.3\r\n\nGenus\t0.9\r\n").unwrap();
    assert_eq!(red.get(Rank::Phylum), Ok(RedValue::new(0.3)));
    assert_eq!(red.get(Rank::Genus), Ok(RedValue::new(0.9)));
}

#[test]
fn full_precision_is_kept() {
    let red = Red::from_tsv("rank\tmedian\nPhylum\t0.23103292252099542\n").unwrap();
    assert_eq!(red.get(Rank::Phylum), Ok(RedValue::new(0.231_032_922_520_995_42)));
}

// Expected value from Python: `RED({...}).to_str()` is
// `d=0.000, p=0.300, c=0.500, o=0.700, f=0.800, g=0.900, s=1.000`.
#[test]
fn to_str_is_in_rank_order_with_three_decimals() {
    assert_eq!(
        red().to_str(),
        "d=0.000, p=0.300, c=0.500, o=0.700, f=0.800, g=0.900, s=1.000"
    );

    // The order of the rows in the input does not matter.
    let red = Red::from_tsv("rank\tmedian\nGenus\t0.9\nPhylum\t0.3\n").unwrap();
    assert_eq!(red.to_str(), "d=0.000, p=0.300, g=0.900, s=1.000");
}

#[test]
fn to_str_rounds_like_python() {
    // Python `'%.3f' % x`: 0.0625 -> 0.062, 0.0635 -> 0.064, 0.9995 -> 1.000.
    let red = Red::from_tsv("rank\tmedian\nPhylum\t0.0625\nClass\t0.0635\nOrder\t0.9995\n").unwrap();
    assert_eq!(red.to_str(), "d=0.000, p=0.062, c=0.064, o=1.000, s=1.000");
}

#[test]
fn from_file_reads_the_tsv() {
    let path = std::env::temp_dir().join(format!("curatomatic_red_{}.tsv", std::process::id()));
    fs::write(&path, TSV).unwrap();
    let from_file = Red::from_file(&path);
    fs::remove_file(&path).unwrap();
    assert_eq!(from_file.unwrap(), red());
}

#[test]
fn from_file_reports_the_path_of_a_bad_file() {
    let path = std::env::temp_dir().join(format!("curatomatic_red_bad_{}.tsv", std::process::id()));
    fs::write(&path, "rank\tmean\n").unwrap();
    let from_file = Red::from_file(&path);
    fs::remove_file(&path).unwrap();
    let message = red_error(from_file);
    assert!(message.contains(&path.display().to_string()), "{message}");
    assert!(message.contains("missing column: median"), "{message}");
}

#[test]
fn from_file_with_a_missing_file_is_an_io_error() {
    assert!(matches!(
        Red::from_file(Path::new("does/not/exist.tsv")),
        Err(CuratomaticError::Io(_))
    ));
}

#[test]
fn from_file_reads_the_r237_file() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("test_data/r237/median_red_rank_distribution.tsv");
    let red = Red::from_file(&path).unwrap();
    let expected = [0.0, 0.3095, 0.4428, 0.5962, 0.7519, 0.9189, 1.0];
    for (rank, value) in Rank::iter().zip(expected) {
        assert_eq!(red.get(rank), Ok(RedValue::new(value)), "{rank:?}");
    }
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
