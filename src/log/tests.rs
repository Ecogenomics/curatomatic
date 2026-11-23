use super::*;

#[test]
fn level_filter_mapping() {
    assert_eq!(LogLevel::Debug.to_filter(), LevelFilter::Debug);
    assert_eq!(LogLevel::Info.to_filter(), LevelFilter::Info);
    assert_eq!(LogLevel::Warning.to_filter(), LevelFilter::Warn);
    assert_eq!(LogLevel::Error.to_filter(), LevelFilter::Error);
}

#[test]
fn critical_maps_to_error_as_the_log_crate_has_no_higher_level() {
    assert_eq!(LogLevel::Critical.to_filter(), LevelFilter::Error);
}

// The Python format is `[%(asctime)s] [%(levelname)s] - %(message)s`, where `WARNING`, `INFO` and
// `CRITICAL` are rewritten to ` WARN`, ` INFO` and ` CRIT` (`curatomatic/util/log.py`).
#[test]
fn line_format_right_aligns_the_level_to_five_characters() {
    let ts = "2026-11-09 12:30:45";
    let cases = [
        (Level::Debug, "[2026-11-09 12:30:45] [DEBUG] - hello\n"),
        (Level::Info, "[2026-11-09 12:30:45] [ INFO] - hello\n"),
        (Level::Warn, "[2026-11-09 12:30:45] [ WARN] - hello\n"),
        (Level::Error, "[2026-11-09 12:30:45] [ERROR] - hello\n"),
    ];
    for (level, expected) in cases {
        assert_eq!(format_line(ts, level, &"hello"), expected);
    }
}

#[test]
fn line_format_formats_the_message_arguments() {
    assert_eq!(
        format_line("t", Level::Info, &format_args!("Found {} nodes", 3)),
        "[t] [ INFO] - Found 3 nodes\n"
    );
}

#[test]
fn timestamp_format_is_twenty_two_characters_with_the_bracket_and_space() {
    let time = chrono::NaiveDate::from_ymd_opt(2026, 11, 9)
        .unwrap()
        .and_hms_opt(1, 2, 3)
        .unwrap();
    let timestamp = time.format(TIMESTAMP_FORMAT).to_string();
    assert_eq!(timestamp, "2026-11-09 01:02:03");
    // The parity check strips the first 22 characters: `[<timestamp>] `.
    let line = format_line(&timestamp, Level::Info, &"x");
    assert_eq!(
        line.chars().take(22).collect::<String>(),
        "[2026-11-09 01:02:03] "
    );
    assert_eq!(line.chars().skip(22).collect::<String>(), "[ INFO] - x\n");
}
