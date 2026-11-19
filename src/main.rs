mod bootstrap;
mod error;
mod genome;
mod log;
mod metadata;
mod newick;
mod node_label;
mod rank;
mod red;
mod report;
mod taxon;
mod taxonomy;
mod taxonomy_file;
mod tree;

use std::collections::BTreeSet;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use ::log::{Level, debug, error, info, log, warn};
use clap::Parser;
use strum::IntoEnumIterator;
use thousands::Separable;

use crate::bootstrap::BootstrapValue;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::log::LogLevel;
use crate::metadata::{COL_ACCESSION, Metadata};
use crate::rank::Rank;
use crate::red::Red;
use crate::report::{ReportInput, build_report, write_report};
use crate::taxon::Taxon;
use crate::taxonomy_file::TaxonomyFile;
use crate::tree::{Tree, fill_in_missing_ranks};

/// Semi-automatic curation of decorated phylogenetic trees.
#[derive(Debug, Parser)]
#[command(version, about)]
struct Cli {
    /// Path to the RED decorated and scaled tree.
    tree: PathBuf,

    /// Path to the median RED of each rank: a TSV with `rank` and `median` columns
    /// (e.g. `median_red_rank_distribution.tsv`).
    #[arg(value_name = "MEDIAN_RED")]
    red_dict: PathBuf,

    /// Output directory to write to.
    out_dir: PathBuf,

    /// Output logging messages equal, or above this level.
    #[arg(long, value_enum, default_value_t = LogLevel::Info)]
    log: LogLevel,

    /// Minimum bootstrap value to consider for decoration.
    #[arg(long, default_value_t = BootstrapValue::new(95.0))]
    min_bs: BootstrapValue,

    /// Tab separated file with no header giving the taxonomy of each genome (e.g. the
    /// `decorate_classifications.tsv` of gDecorate). It names the ranks the tree can't hold, and
    /// is the prior taxonomy for downward refinement.
    #[arg(long, value_name = "PATH")]
    taxonomy: Option<PathBuf>,

    /// Tab separated file with a header (e.g. the GTDB metadata file) giving the name used for
    /// novel taxa of each genome, in the `--novel-accession` and `--novel-name-column` columns.
    #[arg(long, value_name = "PATH", requires = "novel_name_column")]
    novel_name_file: Option<PathBuf>,

    /// Column of `--novel-name-file` holding the genome IDs (e.g. `db_name`).
    #[arg(long, value_name = "NAME", default_value = COL_ACCESSION, requires = "novel_name_file")]
    novel_accession: String,

    /// Column of `--novel-name-file` holding the names (e.g. `ncbi_wgs_formatted`).
    #[arg(long, value_name = "NAME", requires = "novel_name_file")]
    novel_name_column: Option<String>,
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
    debug!("median_red: {}", cli.red_dict.display());
    debug!("out_dir: {}", cli.out_dir.display());
    debug!("min_bs: {}", cli.min_bs);
    debug!("taxonomy: {:?}", cli.taxonomy);
    debug!("novel_name_file: {:?}", cli.novel_name_file);
    debug!("novel_accession: {}", cli.novel_accession);
    debug!("novel_name_column: {:?}", cli.novel_name_column);

    if let Err(e) = curate(cli) {
        error!("{e}");
        return Ok(ExitCode::FAILURE);
    }
    Ok(ExitCode::SUCCESS)
}

/// Log lines that report how many tree genomes have a prior taxonomy, warning when downward
/// refinement can't run because none do.
fn taxonomy_coverage_messages(metadata: &Metadata, tree_n_nodes: usize) -> Vec<(Level, String)> {
    let mut out = Vec::new();
    if !metadata.has_prior_taxonomy {
        out.push((
            Level::Warn,
            "No prior taxonomy was given, downward refinement is disabled.".to_owned(),
        ));
    }
    let n_taxonomy = metadata.n_with_taxonomy();
    out.push((
        Level::Info,
        format!(
            "{} of {} tree genomes have a prior taxonomy.",
            n_taxonomy.separate_with_commas(),
            tree_n_nodes.separate_with_commas()
        ),
    ));
    if metadata.has_prior_taxonomy && n_taxonomy == 0 {
        out.push((
            Level::Warn,
            "No tree genomes have a prior taxonomy, downward refinement is disabled.".to_owned(),
        ));
    }
    out
}

/// Read the taxonomy file, which must hold every genome on the tree.
fn read_taxonomy_file(path: &Path, tree: &Tree) -> CuratomaticResult<TaxonomyFile> {
    info!("Reading taxonomy file");
    let file = TaxonomyFile::from_file(path)?;
    let (missing, extra) = file.compare(&tree.gids);
    if let Some(first) = missing.first() {
        return Err(CuratomaticError::MissingFromTaxonomyFile {
            n: missing.len(),
            first: first.clone(),
        });
    }
    if !extra.is_empty() {
        warn!(
            "Found {} genomes in the taxonomy file that are not in the tree.",
            extra.len().separate_with_commas()
        );
    }
    Ok(file)
}

/// Give the leaves the ranks that the tree can't hold, from the taxonomy file.
fn assign_leaf_taxonomy(tree: &mut Tree, file: &TaxonomyFile) -> CuratomaticResult<()> {
    info!("Assigning ranks from the taxonomy file to the leaves");
    let taxa_before = tree.extract_all_taxa_from_tree()?;
    let assignment = tree.assign_leaf_taxonomy(file)?;
    if !tree.extract_all_taxa_from_tree()?.is_superset(&taxa_before) {
        return Err(CuratomaticError::Invariant(
            "assigning the taxonomy file to the leaves removed a taxon from the tree".to_owned(),
        ));
    }
    let n_taxa = assignment
        .n_taxa
        .iter()
        .map(|(rank, n)| format!("{}__={}", rank.prefix(), n.separate_with_commas()))
        .collect::<Vec<_>>()
        .join(", ");
    info!(
        "Assigned taxa to {} leaves ({n_taxa}).",
        assignment.n_leaves.separate_with_commas()
    );
    if !assignment.skipped.is_empty() {
        warn!(
            "Skipped {} taxa of the taxonomy file that are already on the tree, or that more than \
             one leaf would get.",
            assignment.skipped.len().separate_with_commas()
        );
        for taxon in &assignment.skipped {
            debug!("Skipped: {taxon}");
        }
    }
    Ok(())
}

/// Read the novel names and the prior taxonomy (from the taxonomy file), subset to the genomes on
/// the tree, and every taxon of the prior taxonomy for the report. Each is empty if its file is
/// not given.
fn read_metadata(
    cli: &Cli,
    tree: &Tree,
    tree_n_nodes: usize,
    taxonomy: Option<&TaxonomyFile>,
) -> CuratomaticResult<(Metadata, BTreeSet<Taxon>)> {
    let mut metadata = Metadata::default();
    if let (Some(path), Some(column)) = (&cli.novel_name_file, &cli.novel_name_column) {
        info!("Reading novel names");
        metadata = Metadata::from_file(path, &cli.novel_accession, column)?;
        let n_meta_rows_before = metadata.rows.len();
        metadata.subset_to_ids(&tree.gids);
        let n_meta_mismatch = tree_n_nodes.saturating_sub(metadata.rows.len());
        info!(
            "Found {} records, matched {} IDs to the tree.",
            n_meta_rows_before.separate_with_commas(),
            metadata.rows.len().separate_with_commas()
        );
        if n_meta_mismatch > 0 {
            warn!(
                "Found {} IDs in the tree that do not have a novel name entry.",
                n_meta_mismatch.separate_with_commas()
            );
        }
    }

    for (gid, gid_taxonomy) in taxonomy.into_iter().flat_map(TaxonomyFile::iter) {
        metadata.set_taxonomy(gid.clone(), gid_taxonomy.clone());
    }
    let prior_taxa = metadata.all_taxa();
    metadata.subset_to_ids(&tree.gids);
    for (level, message) in taxonomy_coverage_messages(&metadata, tree_n_nodes) {
        log!(level, "{message}");
    }
    Ok((metadata, prior_taxa))
}

/// Check that the taxa on the tree are those it started with, plus the novel taxa.
fn check_taxa_unchanged(
    before: &BTreeSet<Taxon>,
    tree: &Tree,
    novel: &BTreeSet<Taxon>,
) -> CuratomaticResult<()> {
    let after = tree.extract_all_taxa_from_tree()?;
    let lost: Vec<String> = before
        .iter()
        .filter(|taxon| !after.contains(*taxon) || novel.contains(*taxon))
        .map(ToString::to_string)
        .collect();
    let gained: Vec<String> = after
        .iter()
        .filter(|taxon| !before.contains(*taxon) && !novel.contains(*taxon))
        .map(ToString::to_string)
        .collect();
    if lost.is_empty() && gained.is_empty() {
        return Ok(());
    }
    Err(CuratomaticError::Invariant(format!(
        "the taxa on the tree changed by more than the novel taxa that were created \
         (lost or reused by a novel taxon: [{}], gained: [{}])",
        lost.join(", "),
        gained.join(", ")
    )))
}

fn curate(cli: &Cli) -> CuratomaticResult<()> {
    let stem = cli.tree.file_stem().unwrap_or_default().to_string_lossy();
    let path_tree_out = cli.out_dir.join(format!("{stem}.curatomatic.tree"));
    let path_report_out = cli.out_dir.join(format!("{stem}.curatomatic.report.tsv"));
    let path_taxonomy_out = cli.out_dir.join(format!("{stem}.curatomatic.taxonomy.tsv"));

    info!("Reading tree");
    let mut tree = Tree::from_file(&cli.tree)?;

    let tree_n_nodes = tree.n_leaves();
    info!(
        "Found: {} leaf nodes in the tree.",
        tree_n_nodes.separate_with_commas()
    );

    // This runs first, before anything reads the taxa on the tree.
    let taxonomy = match &cli.taxonomy {
        Some(path) => {
            let file = read_taxonomy_file(path, &tree)?;
            assign_leaf_taxonomy(&mut tree, &file)?;
            Some(file)
        }
        None => None,
    };

    let tree_taxa_before = tree.extract_all_taxa_from_tree()?;
    let taxon_nodes_before = tree.taxon_nodes();

    let (metadata, prior_taxa) = read_metadata(cli, &tree, tree_n_nodes, taxonomy.as_ref())?;
    drop(taxonomy);
    // Read the median RED of each rank.
    info!("Reading median RED values");
    let red_dict = Red::from_file(&cli.red_dict)?;
    debug!("{}", red_dict.to_str());

    let mut moves = Vec::new();
    for rank in Rank::iter() {
        if matches!(rank, Rank::Domain | Rank::Species) {
            continue;
        }
        info!("Refining existing rank: {}", rank.prefix());
        let refinement = tree.refine_existing_labels(rank, cli.min_bs, &red_dict, &metadata)?;
        moves.extend(tree.process_refinement(&refinement, rank, cli.min_bs, &red_dict)?);
        tree.reload()?;
    }

    info!("Identifying leaf nodes missing ranks");
    let leaf_to_ranks_missing = tree.get_leaf_nodes_missing_rank()?;
    info!(
        "Found {} genomes missing one or more ranks.",
        leaf_to_ranks_missing.len().separate_with_commas()
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
    info!(
        "Created {} novel ranks.",
        new_taxa_created.len().separate_with_commas()
    );

    check_taxa_unchanged(&tree_taxa_before, &tree, &new_taxa_created)?;

    info!("Building report");
    tree.reload()?;
    let report = build_report(&ReportInput {
        tree: &tree,
        taxa_before: &taxon_nodes_before,
        moves: &moves,
        novel: &new_taxa_created,
        metadata: &metadata,
        prior_taxa: &prior_taxa,
        min_bs: cli.min_bs,
        red_dict: &red_dict,
    })?;

    info!("Writing taxonomy to disk.");
    let leaf_taxonomies = tree.leaf_taxonomies()?;
    taxonomy_file::write(
        &path_taxonomy_out,
        leaf_taxonomies
            .iter()
            .map(|(gid, taxa)| (gid, taxa.as_slice())),
    )?;

    info!("Setting the output labels (leaves by genome ID, no RED values)");
    tree.set_output_labels()?;

    info!("Writing tree to disk.");
    tree.to_file(&path_tree_out)?;

    info!(
        "Writing report to disk ({} rows).",
        report.len().separate_with_commas()
    );
    write_report(&path_report_out, &report)?;

    info!("Done.");
    Ok(())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    run(&cli).unwrap_or_else(|e| {
        // Ignore a failure to write to stderr, there is nothing else to report it to.
        drop(writeln!(io::stderr(), "ERROR: {e}"));
        ExitCode::FAILURE
    })
}

#[cfg(test)]
#[path = "main/tests.rs"]
mod tests;
