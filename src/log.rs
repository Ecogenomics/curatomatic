use std::fmt;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use std::sync::Mutex;

use ::log::{Level, LevelFilter, Log, Metadata, Record};
use clap::ValueEnum;

/// The verbosity of the logger.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub(crate) enum LogLevel {
    Debug,
    Info,
    Warning,
    Error,
    Critical,
}

impl LogLevel {
    /// The `log` crate has no level above error, so `critical` maps to `error`.
    fn to_filter(self) -> LevelFilter {
        match self {
            LogLevel::Debug => LevelFilter::Debug,
            LogLevel::Info => LevelFilter::Info,
            LogLevel::Warning => LevelFilter::Warn,
            LogLevel::Error | LogLevel::Critical => LevelFilter::Error,
        }
    }
}

/// Name of the log file written within the output directory.
const LOG_FILE_NAME: &str = "curatomatic.log";

/// The format of the timestamp at the start of each line.
const TIMESTAMP_FORMAT: &str = "%Y-%m-%d %H:%M:%S";

/// Format a log line as `[timestamp] [LEVEL] - message`, with the level right-aligned to 5
/// characters (e.g. ` INFO`, ` WARN`, `ERROR`) like the Python implementation.
fn format_line(timestamp: &str, level: Level, message: &impl fmt::Display) -> String {
    format!("[{timestamp}] [{level:>5}] - {message}\n")
}

/// Logs to both stderr and a file, formatted as `[time] [LEVEL] - message`.
struct Logger {
    level: LevelFilter,
    file: Mutex<File>,
}

impl Log for Logger {
    fn enabled(&self, metadata: &Metadata<'_>) -> bool {
        metadata.level() <= self.level
    }

    fn log(&self, record: &Record<'_>) {
        if !self.enabled(record.metadata()) {
            return;
        }

        let timestamp = chrono::Local::now().format(TIMESTAMP_FORMAT).to_string();
        let line = format_line(&timestamp, record.level(), record.args());

        // A logger has nowhere to report its own write failures, so they are dropped.
        drop(io::stderr().lock().write_all(line.as_bytes()));
        if let Ok(mut file) = self.file.lock() {
            drop(file.write_all(line.as_bytes()));
        }
    }

    fn flush(&self) {
        drop(io::stderr().flush());
        if let Ok(mut file) = self.file.lock() {
            drop(file.flush());
        }
    }
}

/// Install the global logger, writing to stderr and `<out_dir>/curatomatic.log`.
///
/// The output directory must already exist.
pub(crate) fn init(out_dir: &Path, level: LogLevel) -> io::Result<()> {
    let level = level.to_filter();
    let file = File::create(out_dir.join(LOG_FILE_NAME))?;
    log::set_boxed_logger(Box::new(Logger {
        level,
        file: Mutex::new(file),
    }))
    .map_err(io::Error::other)?;
    log::set_max_level(level);
    Ok(())
}

#[cfg(test)]
mod tests;
