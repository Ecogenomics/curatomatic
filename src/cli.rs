use std::path::PathBuf;

use clap::Parser;

use crate::bootstrap::BootstrapValue;
use crate::log::LogLevel;
use crate::metadata::COL_ACCESSION;

/// Semi-automatic curation of decorated phylogenetic trees.
#[derive(Debug, Parser)]
#[command(version, about)]
pub(crate) struct Cli {
    /// Path to the RED decorated and scaled tree.
    pub(crate) tree: PathBuf,

    /// Path to the median RED of each rank: a TSV with `rank` and `median` columns
    /// (e.g. `median_red_rank_distribution.tsv`).
    #[arg(value_name = "MEDIAN_RED")]
    pub(crate) red_dict: PathBuf,

    /// Output directory to write to.
    pub(crate) out_dir: PathBuf,

    /// Output logging messages equal, or above this level.
    #[arg(long, value_enum, default_value_t = LogLevel::Info)]
    pub(crate) log: LogLevel,

    /// Minimum bootstrap value to consider for decoration.
    #[arg(long, default_value_t = BootstrapValue::new(95.0))]
    pub(crate) min_bs: BootstrapValue,

    /// Tab separated file with no header giving the taxonomy of each genome (e.g. the
    /// `decorate_classifications.tsv` of gDecorate). It names the ranks the tree can't hold, and
    /// is the prior taxonomy for downward refinement. Required: it gives every leaf its species.
    #[arg(long, value_name = "PATH")]
    pub(crate) taxonomy: PathBuf,

    /// Tab separated file with a header (e.g. the GTDB metadata file) giving the name used for
    /// novel taxa of each genome, in the `--novel-accession` and `--novel-name-column` columns.
    #[arg(long, value_name = "PATH", requires = "novel_name_column")]
    pub(crate) novel_name_file: Option<PathBuf>,

    /// Column of `--novel-name-file` holding the genome IDs (e.g. `db_name`).
    #[arg(long, value_name = "NAME", default_value = COL_ACCESSION, requires = "novel_name_file")]
    pub(crate) novel_accession: String,

    /// Column of `--novel-name-file` holding the names (e.g. `ncbi_wgs_formatted`).
    #[arg(long, value_name = "NAME", requires = "novel_name_file")]
    pub(crate) novel_name_column: Option<String>,

    /// Also write the inputs as curatomatic loaded them, before any curation: the tree with a
    /// dummy node above each leaf holding the leaf's taxa (`<stem>.debug.zombie.tree`), and the
    /// full taxonomy of each genome (`<stem>.debug.taxonomy.tsv`).
    #[arg(long)]
    pub(crate) debug: bool,
}

#[cfg(test)]
mod tests;
