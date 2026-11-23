mod bootstrap;
mod cli;
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
use std::path::Path;
use std::process::ExitCode;

use ::log::{Level, debug, error, info, log, warn};
use clap::Parser;
use strum::IntoEnumIterator;
use thousands::Separable;

use crate::cli::Cli;
use crate::error::{CuratomaticError, CuratomaticResult};
use crate::log::LogLevel;
use crate::metadata::Metadata;
use crate::rank::Rank;
use crate::red::Red;
use crate::report::{ReportInput, build_report, write_report};
use crate::taxon::{Taxon, TaxonCounts};
use crate::taxonomy_file::TaxonomyFile;
use crate::tree::{Tree, fill_in_missing_ranks};

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
    debug!("taxonomy: {}", cli.taxonomy.display());
    debug!("prior_taxonomy: {:?}", cli.prior_taxonomy);
    debug!("novel_name_file: {:?}", cli.novel_name_file);
    debug!("novel_accession: {}", cli.novel_accession);
    debug!("novel_name_column: {:?}", cli.novel_name_column);
    debug!("debug: {}", cli.debug);

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

/// Give the leaves the taxa that the tree can't hold, from the taxonomy file.
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
    Ok(())
}

/// Read the novel names and the prior taxonomy (`--prior-taxonomy`), subset to the genomes on the
/// tree, and every taxon of the prior taxonomy for the report. The novel names are empty if their
/// file is not given, and no genome has a prior taxonomy if `prior` is `None`.
///
/// The prior taxonomy is never read from `--taxonomy`: that is the decoration of the same tree,
/// so every genome below a taxon would have it, and a taxon could never be moved down.
fn read_metadata(
    cli: &Cli,
    tree: &Tree,
    tree_n_nodes: usize,
    prior: Option<&TaxonomyFile>,
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

    for (gid, gid_taxonomy) in prior.iter().flat_map(|prior| prior.iter()) {
        metadata.set_taxonomy(gid.clone(), gid_taxonomy.clone());
    }
    let prior_taxa = metadata.all_taxa();
    metadata.subset_to_ids(&tree.gids);
    for (level, message) in taxonomy_coverage_messages(&metadata, tree_n_nodes) {
        log!(level, "{message}");
    }
    Ok((metadata, prior_taxa))
}

/// Check that curatomatic only moved or created taxa (see the rules in `AGENTS.md`).
///
/// `before` counts the named taxa on the tree before refinement, and `after` the named taxa at the
/// end of the run. Every taxon of `before` must still be on the tree, on as many nodes as before
/// (it was moved, not removed or copied), with its text unchanged. Every other taxon on the tree
/// must be one of the `novel` taxa that were created, and a novel taxon must not reuse a name that
/// was already on the tree. Placeholders (e.g. `g__`) are not counted, they can be refined freely.
fn check_taxa_unchanged(
    before: &TaxonCounts,
    after: &TaxonCounts,
    novel: &BTreeSet<Taxon>,
) -> CuratomaticResult<()> {
    let join = |taxa: Vec<String>| taxa.join(", ");
    let removed: Vec<String> = before
        .iter()
        .filter(|&(taxon, _)| after.get(taxon) == 0)
        .map(|(taxon, _)| taxon.to_string())
        .collect();
    let count_changed: Vec<String> = before
        .iter()
        .filter(|&(taxon, n)| after.get(taxon) != 0 && after.get(taxon) != n)
        .map(|(taxon, n)| format!("{taxon} ({n} -> {})", after.get(taxon)))
        .collect();
    let reused: Vec<String> = novel
        .iter()
        .filter(|taxon| before.get(taxon) != 0)
        .map(ToString::to_string)
        .collect();
    let gained: Vec<String> = after
        .iter()
        .filter(|&(taxon, _)| before.get(taxon) == 0 && !novel.contains(taxon))
        .map(|(taxon, _)| taxon.to_string())
        .collect();

    let mut problems = Vec::new();
    if !removed.is_empty() {
        problems.push(format!("removed: [{}]", join(removed)));
    }
    if !count_changed.is_empty() {
        problems.push(format!(
            "on a different number of nodes: [{}]",
            join(count_changed)
        ));
    }
    if !reused.is_empty() {
        problems.push(format!("reused by a novel taxon: [{}]", join(reused)));
    }
    if !gained.is_empty() {
        problems.push(format!(
            "gained, but not created as a novel taxon: [{}]",
            join(gained)
        ));
    }
    if problems.is_empty() {
        return Ok(());
    }
    Err(CuratomaticError::TaxaChanged(problems.join("; ")))
}

/// Write the inputs as they were loaded (`--debug`): the tree with a dummy node above each leaf,
/// and the full taxonomy of each genome. Removing the dummy nodes of an old style tree and
/// loading it with its taxonomy file must give these back unchanged.
fn write_debug_inputs(tree: &Tree, out_dir: &Path, stem: &str) -> CuratomaticResult<()> {
    let path_tree = out_dir.join(format!("{stem}.debug.zombie.tree"));
    let path_taxonomy = out_dir.join(format!("{stem}.debug.taxonomy.tsv"));
    info!("Writing the loaded inputs with dummy nodes (--debug).");
    std::fs::write(&path_tree, tree.to_zombie_newick())
        .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path_tree.display())))?;
    let taxonomies = tree.full_taxonomies()?;
    taxonomy_file::write(
        &path_taxonomy,
        taxonomies.iter().map(|(gid, taxa)| (gid, taxa.as_slice())),
    )
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
    let taxonomy = read_taxonomy_file(&cli.taxonomy, &tree)?;
    assign_leaf_taxonomy(&mut tree, &taxonomy)?;
    // The inputs are gospel: a leaf without its species or domain is an error, not something to
    // curate.
    tree.check_leaves_have_species()?;
    tree.check_leaves_have_domain()?;

    if cli.debug {
        write_debug_inputs(&tree, &cli.out_dir, &stem)?;
    }

    // Every change to the taxa from here on must be a move or a novel taxon.
    let tree_taxa_before = tree.taxon_counts()?;

    drop(taxonomy);
    let prior = match &cli.prior_taxonomy {
        Some(path) => {
            info!("Reading prior taxonomy");
            Some(TaxonomyFile::from_file(path)?)
        }
        None => None,
    };
    let (metadata, prior_taxa) = read_metadata(cli, &tree, tree_n_nodes, prior.as_ref())?;
    drop(prior);
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

    // Refinement and novel taxa are the only steps that change the taxa, so this is checked
    // before anything is written: a run that breaks the rule fails without any output.
    info!("Checking that taxa were only moved or created");
    check_taxa_unchanged(&tree_taxa_before, &tree.taxon_counts()?, &new_taxa_created)?;
    info!("Checking that every leaf holds its species and has a domain");
    tree.check_leaves_have_species()?;
    tree.check_leaves_have_domain()?;

    info!("Building report");
    tree.reload()?;
    let report = build_report(&ReportInput {
        tree: &tree,
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

    info!("Removing RED labels");
    tree.remove_red_labels()?;

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
