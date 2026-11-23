//! Opening the input files, which may be gzipped.
//!
//! An input whose name ends in `.gz` is decompressed as it is read, and any other input is read
//! as plain text. This is decided by the name only, never by the content. Every input is opened
//! with [`open_input`] (or [`read_input_to_string`]), never with `File::open` or
//! `fs::read_to_string` directly.

use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::{CuratomaticError, CuratomaticResult};

/// Whether `path` is read gzipped: its name ends in `.gz`.
fn is_gzipped(path: &Path) -> bool {
    path.extension() == Some(OsStr::new("gz"))
}

/// The error of an input that could not be read, naming its path.
fn io_error(path: &Path, error: &std::io::Error) -> CuratomaticError {
    CuratomaticError::Io(format!("{}: {error}", path.display()))
}

/// Open an input file, decompressing it if its name ends in `.gz` (see the module docs).
///
/// A gzipped file is read with [`MultiGzDecoder`], so every member of a bgzip or concatenated
/// gzip file is read, not only the first. A file that is not valid gzip, or is cut short, gives an
/// error when it is read, which the caller reports as a [`CuratomaticError::Io`] naming the path.
pub(crate) fn open_input(path: &Path) -> CuratomaticResult<Box<dyn BufRead>> {
    #[expect(
        clippy::disallowed_methods,
        reason = "the one place an input is opened, see the module docs"
    )]
    let file = File::open(path).map_err(|e| io_error(path, &e))?;
    if is_gzipped(path) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(
            BufReader::new(file),
        ))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Read the whole of an input file into a string, see [`open_input`].
pub(crate) fn read_input_to_string(path: &Path) -> CuratomaticResult<String> {
    let mut text = String::new();
    open_input(path)?
        .read_to_string(&mut text)
        .map_err(|e| io_error(path, &e))?;
    Ok(text)
}

/// The stem of an input, used to name the outputs: its file name without a trailing `.gz`, and
/// then without its last extension. `x.nwk.gz` and `x.nwk` both give `x`.
pub(crate) fn input_stem(path: &Path) -> String {
    let path = if is_gzipped(path) {
        path.with_extension("")
    } else {
        path.to_path_buf()
    };
    path.file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned()
}

#[cfg(test)]
mod tests;
