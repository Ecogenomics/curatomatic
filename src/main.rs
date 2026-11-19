mod bootstrap;
mod error;
mod genome;
mod log;
mod metadata;
mod newick;
mod node_label;
mod rank;
mod red;
mod taxon;
mod taxonomy;
mod tree;

use std::collections::BTreeSet;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use ::log::{Level, debug, error, info, log, warn};
use clap::Parser;
use strum::IntoEnumIterator;

use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::log::LogLevel;
use crate::metadata::Metadata;
use crate::rank::Rank;
use crate::red::Red;
use crate::tree::{Tree, fill_in_missing_ranks};

/// Semi-automatic curation of decorated phylogenetic trees.
#[derive(Debug, Parser)]
#[command(version, about)]
struct Cli {
    /// Path to the RED decorated and scaled tree.
    tree: PathBuf,

    /// Path to the RED dictionary.
    red_dict: PathBuf,

    /// Output directory to write to.
    out_dir: PathBuf,

    /// Output logging messages equal, or above this level.
    #[arg(long, value_enum, default_value_t = LogLevel::Info)]
    log: LogLevel,

    /// Minimum bootstrap value to consider for decoration.
    #[arg(long, default_value_t = BootstrapValue::new(95.0))]
    min_bs: BootstrapValue,

    /// File containing unique identifiers for novel rank creation.
    #[arg(long)]
    meta: Option<PathBuf>,
}

fn run(cli: &Cli) -> io::Result<ExitCode> {
    if cli.log != LogLevel::Debug && cli.out_dir.exists() {
        writeln!(
            io::stderr(),
            "ERROR: Output directory already exists: {}",
            cli.out_dir.display()
        )?;
        return Ok(ExitCode::FAILURE);
    }

    std::fs::create_dir_all(&cli.out_dir)?;
    log::init(&cli.out_dir, cli.log)?;

    info!(
        "curatomatic v{} {}",
        env!("CARGO_PKG_VERSION"),
        std::env::args().skip(1).collect::<Vec<_>>().join(" ")
    );
    debug!("tree: {}", cli.tree.display());
    debug!("red_dict: {}", cli.red_dict.display());
    debug!("out_dir: {}", cli.out_dir.display());
    debug!("min_bs: {}", cli.min_bs);
    debug!("meta: {:?}", cli.meta);

    if let Err(e) = curate(cli) {
        error!("{e}");
        return Ok(ExitCode::FAILURE);
    }
    Ok(ExitCode::SUCCESS)
}

/// Format a count with thousands separators, e.g. `1,234`.
fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::new();
    for (i, digit) in digits.chars().enumerate() {
        if i > 0 && digits.len().saturating_sub(i) % 3 == 0 {
            out.push(',');
        }
        out.push(digit);
    }
    out
}

/// Log lines that report how many tree genomes have a prior taxonomy, warning when downward
/// refinement can't run because none do.
fn taxonomy_coverage_messages(metadata: &Metadata, tree_n_nodes: usize) -> Vec<(Level, String)> {
    let mut out = Vec::new();
    if !metadata.has_taxonomy_column {
        out.push((
            Level::Warn,
            "The metadata file has no gtdb_taxonomy column, downward refinement is disabled."
                .to_owned(),
        ));
    }
    let n_taxonomy = metadata.n_with_taxonomy();
    out.push((
        Level::Info,
        format!(
            "{} of {} tree genomes have a gtdb_taxonomy.",
            thousands(n_taxonomy),
            thousands(tree_n_nodes)
        ),
    ));
    if metadata.has_taxonomy_column && n_taxonomy == 0 {
        out.push((
            Level::Warn,
            "No tree genomes have a gtdb_taxonomy, downward refinement is disabled.".to_owned(),
        ));
    }
    out
}

fn curate(cli: &Cli) -> CuratomaticResult<()> {
    let stem = cli.tree.file_stem().unwrap_or_default().to_string_lossy();
    let path_tree_out = cli.out_dir.join(format!("{stem}.curatomatic.tree"));

    info!("Reading tree");
    let mut tree = Tree::from_file(&cli.tree)?;

    let tree_taxa_before = tree.extract_all_taxa_from_tree()?;

    // Get the number of tree nodes, and zombie nodes (i.e. prefixed with "D-").
    let (tree_n_nodes, tree_n_zombies) = tree.n_nodes_and_zombies()?;
    info!(
        "Found: {} leaf nodes, and {} zombies in the tree.",
        thousands(tree_n_nodes),
        thousands(tree_n_zombies)
    );

    // Read the metadata file.
    let metadata = if let Some(meta) = &cli.meta {
        info!("Reading metadata");
        let mut metadata = Metadata::from_file(meta)?;
        let n_meta_rows_before = metadata.rows.len();
        metadata.subset_to_ids(&tree.gids);
        let n_meta_mismatch = tree_n_nodes.saturating_sub(metadata.rows.len());
        info!(
            "Found {} records, matched {} IDs to the tree.",
            thousands(n_meta_rows_before),
            thousands(metadata.rows.len())
        );
        if n_meta_mismatch > 0 {
            warn!(
                "Found {} IDs in the tree that do not have a metadata entry.",
                thousands(n_meta_mismatch)
            );
        }
        for (level, message) in taxonomy_coverage_messages(&metadata, tree_n_nodes) {
            log!(level, "{message}");
        }
        metadata
    } else {
        Metadata::default()
    };

    // Read the RED dictionary.
    info!("Reading RED dictionary");
    let red_dict = Red::from_file(&cli.red_dict)?;
    debug!("{}", red_dict.to_str());

    for rank in Rank::iter() {
        if matches!(rank, Rank::Domain | Rank::Species) {
            continue;
        }
        info!("Refining existing rank: {}", rank.prefix());
        let refinement = tree.refine_existing_labels(rank, cli.min_bs, &red_dict, &metadata)?;
        tree.process_refinement(&refinement, rank, cli.min_bs, &red_dict)?;
        tree.reload()?;
        tree.check_no_same_rank_taxa()?;
    }

    info!("Identifying leaf nodes missing ranks");
    let leaf_to_ranks_missing = tree.get_leaf_nodes_missing_rank()?;
    info!(
        "Found {} genomes missing one or more ranks.",
        thousands(leaf_to_ranks_missing.len())
    );

    info!("Creating novel ranks");
    let known_taxa = tree.taxon_labels.clone();
    let new_taxa_created = fill_in_missing_ranks(
        &mut tree,
        &leaf_to_ranks_missing,
        cli.min_bs,
        &red_dict,
        &metadata,
        known_taxa,
    )?;
    info!("Created {} novel ranks.", thousands(new_taxa_created.len()));
    tree.check_no_same_rank_taxa()?;

    let tree_taxa_after = tree.extract_all_taxa_from_tree()?;
    let taxa_without_new: BTreeSet<_> = tree_taxa_after
        .difference(&new_taxa_created)
        .cloned()
        .collect();
    if taxa_without_new != tree_taxa_before {
        return Err(CuratomaticError::Invariant(
            "the taxa on the tree changed by more than the novel taxa that were created".to_owned(),
        ));
    }

    info!("Removing RED labels");
    tree.remove_red_labels()?;

    info!("Writing tree to disk.");
    tree.to_file(&path_tree_out)?;

    info!("Done.");
    Ok(())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(code) => code,
        Err(e) => {
            // Ignore a failure to write to stderr, there is nothing else to report it to.
            drop(writeln!(io::stderr(), "ERROR: {e}"));
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
#[path = "main/tests.rs"]
mod tests;
