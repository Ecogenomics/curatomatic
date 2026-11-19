use std::error::Error;
use std::fmt;

use crate::genome::GenomeId;

/// Errors that can occur within curatomatic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum CuratomaticError {
    /// A taxon was an empty string, so its rank could not be determined.
    EmptyTaxon,
    /// The first character of the taxon is not a known rank prefix.
    UnknownRank(String),
    /// The taxon to remove is not present in the taxonomy.
    TaxonNotFound(String),
    /// A file could not be read.
    Io(String),
    /// A required column is absent from the header of a file.
    MissingColumn { path: String, column: String },
    /// A row has no value for a column that the header declares.
    MissingValue { line: usize, column: String },
    /// The rank (identified by its prefix) has no parent rank, i.e. it is the domain.
    NoParentRank(char),
    /// A Newick tree could not be parsed.
    Newick(String),
    /// A node label does not have the expected `<bootstrap>:<taxonomy>|RED=<value>` format.
    InvalidNodeLabel(String),
    /// A taxon label is found on more than one node of the tree.
    NonUniqueTaxon(String),
    /// The median RED file could not be parsed.
    Red(String),
    /// The median RED file has no value for a rank (identified by its prefix).
    MissingRed(char),
    /// No node was found that a novel taxon for the rank (identified by its prefix) can be placed on.
    NoCandidateNode { rank: char, gid: GenomeId },
    /// An internal expectation was violated.
    Invariant(String),
}

impl fmt::Display for CuratomaticError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CuratomaticError::EmptyTaxon => {
                write!(f, "Empty taxon, unable to determine its rank")
            }
            CuratomaticError::UnknownRank(taxon) => write!(f, "Unknown rank for taxon: {taxon}"),
            CuratomaticError::TaxonNotFound(taxon) => write!(f, "Taxon not found: {taxon}"),
            CuratomaticError::Io(message) => write!(f, "I/O error: {message}"),
            CuratomaticError::MissingColumn { path, column } => {
                write!(f, "Missing required column {column} in: {path}")
            }
            CuratomaticError::MissingValue { line, column } => {
                write!(f, "Line {line} has no value for column: {column}")
            }
            CuratomaticError::NoParentRank(rank) => {
                write!(f, "The {rank}__ rank has no parent rank.")
            }
            CuratomaticError::Newick(message) => write!(f, "Invalid Newick tree: {message}"),
            CuratomaticError::InvalidNodeLabel(label) => write!(f, "Invalid node label: {label}"),
            CuratomaticError::NonUniqueTaxon(taxon) => {
                write!(f, "Taxon label {taxon} is not unique on the tree.")
            }
            CuratomaticError::Red(message) => write!(f, "Invalid median RED file: {message}"),
            CuratomaticError::MissingRed(rank) => {
                write!(f, "The median RED file has no value for rank: {rank}")
            }
            CuratomaticError::NoCandidateNode { rank, gid } => {
                write!(
                    f,
                    "No node can hold a novel {rank}__ taxon for genome: {gid}"
                )
            }
            CuratomaticError::Invariant(message) => write!(f, "Internal error: {message}"),
        }
    }
}

impl Error for CuratomaticError {}

/// A [`Result`] whose error type is [`CuratomaticError`].
pub(crate) type CuratomaticResult<T> = Result<T, CuratomaticError>;
