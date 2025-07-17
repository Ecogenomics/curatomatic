import sys
from pathlib import Path

import typer
from typing_extensions import Annotated

from curatomatic import __version__
from curatomatic.model.log import LogLevel
from curatomatic.model.metadata import Metadata
from curatomatic.model.rank import RANKS
from curatomatic.model.red import RED
from curatomatic.model.report import Report
from curatomatic.model.tree import Tree
from curatomatic.util.log import log_setup
from curatomatic.util.tree import check_for_candidate_sharing, refine_nodes, fill_in_missing_ranks, create_novel_ranks

app = typer.Typer()

TREE_HELP = "Path to the RED decorated and scaled tree."
RED_DICT_HELP = "Path to the RED dictionary."
OUT_DIR_HELP = "Output directory to write to."
LOG_HELP = "Output logging messages equal, or above this level."
MIN_BS_HELP = "Minimum bootstrap value to consider for decoration."
META_HELP = "File containing unique identifiers for novel rank creation."


@app.command()
def run(
        tree: Annotated[Path, typer.Argument(help=TREE_HELP)],
        red_dict: Annotated[Path, typer.Argument(help=RED_DICT_HELP)],
        out_dir: Annotated[Path, typer.Argument(help=OUT_DIR_HELP)],
        log: Annotated[LogLevel, typer.Option(help=LOG_HELP)] = LogLevel.INFO,
        min_bs: Annotated[float, typer.Option(help=MIN_BS_HELP)] = 95.0,
        meta: Annotated[Path, typer.Option(help=META_HELP)] = None,
):
    # Check if the output directory already exists
    if not log is LogLevel.DEBUG and out_dir.exists():
        print(f'ERROR: Output directory already exists: {out_dir}')
        sys.exit(1)

    # Create the output directory and paths
    out_dir.mkdir(parents=True, exist_ok=True)
    path_report_out = out_dir / 'report.tsv'
    path_tree_out = out_dir / f'{tree.stem}.curatomatic.tree'

    # Create the logger
    logger = log_setup(out_dir=out_dir, level=log)
    logger.info(f'curatomatic v{__version__} {" ".join(sys.argv[1:])}')

    # Read the tree
    logger.info(f'Reading tree')
    tree = Tree.from_file(tree)

    tree_taxa_before = tree.extract_all_taxa_from_tree()

    # Get the number of tree nodes, and zombie nodes (i.e. prefixed with "D-").
    tree_n_nodes, tree_n_zombies = tree.n_nodes_and_zombies()
    logger.info(f'Found: {tree_n_nodes:,} leaf nodes, and {tree_n_zombies:,} zombies in the tree.')

    # Read the metadata file
    if meta is not None:
        logger.info(f'Reading metadata')
        metadata = Metadata.from_file(meta)
        n_meta_rows_before = len(metadata.rows)
        metadata.subset_to_ids(set(tree.gid_to_node.keys()))
        n_meta_mismatch = tree_n_nodes - len(metadata.rows)
        logger.info(f'Found {n_meta_rows_before:,} records, matched {len(metadata.rows):,} IDs to the tree.')
        if n_meta_mismatch > 0:
            logger.warning(f'Found {n_meta_mismatch:,} IDs in the tree that do not have a metadata entry.')
    else:
        metadata = Metadata(dict())

    # Read the RED dictionary
    logger.info(f'Reading RED dictionary')
    red_dict = RED.from_file(red_dict)
    logger.debug(red_dict.to_str())

    for rank in RANKS:
        if rank in {'d', 's'}:
            continue
        logger.info(f'Refining existing rank: {rank}')
        cur_rank_refinement = tree.refine_existing_labels(rank, min_bs, red_dict, metadata)
        tree.process_refinement(cur_rank_refinement, rank, min_bs, red_dict, metadata)
        tree.reload()

    logger.info('Identifying leaf nodes missing ranks')
    d_leaf_to_ranks_missing = tree.get_leaf_nodes_missing_rank()
    logger.info(f'Found {len(d_leaf_to_ranks_missing):,} genomes missing one or more ranks.')

    logger.info('Creating novel ranks')
    new_taxa_created = fill_in_missing_ranks(
        d_leaf_to_ranks_missing,
        tree.node_to_descs,
        min_bs,
        red_dict,
        metadata,
        set(tree.taxon_labels_to_node.keys())
    )
    logger.info(f'Created {len(new_taxa_created):,} novel ranks.')

    tree_taxa_after = tree.extract_all_taxa_from_tree()
    assert tree_taxa_after - new_taxa_created == tree_taxa_before

    logger.info('Removing RED labels')
    tree.remove_red_labels()

    # logger.info('Writing report to disk.')
    # report = Report(moved, new_taxa_created)
    # report.write(path_report_out)

    logger.info('Writing tree to disk.')
    tree.to_file(path_tree_out)

    logger.info('Done.')
    sys.exit(0)

    return


    """
    OLD CODE BELOW, IGNORE
    """


    """
    Perform an initial refinement of the tree to check if pre-existing labels
    can be moved to a node with higher bootstrap support.
    
    """
    # tree.find_low_bootstrap_support_nodes_with_taxon_placement(min_bs, red_dict, metadata)



    """
    The tree must have zombies, if not add them.
    
    
    Step 1. Find what nodes can be refined (moving them up, to include new genomes)
            or have better support. These are never removed.
            
    Step 2. Create de novo ranks
    
    """
    poly = tree.get_nested_ranks()
    if len(poly) > 0:
        logger.warning(f'Found {len(poly):,} nested ranks, these taxa will be ignored.')

    d_taxon_to_candidates = tree.find_taxa_for_refinement(poly=poly,min_bootstrap=min_bs, red_dict=red_dict)
    logger.info(f'Found {len(d_taxon_to_candidates):,} existing taxa that can be moved higher.')

    # Check if there are any taxa that share a candidate node
    shared = check_for_candidate_sharing(d_taxon_to_candidates)
    if len(shared) > 0:
        logger.warning(f'There are {len(shared):,} taxa that share a candidate node.')
        raise NotImplemented("!")

    # Move those candidate nodes up
    moved = refine_nodes(d_taxon_to_candidates, tree.taxon_labels_to_node)
    moved = set(moved)
    logger.info(f'Moved {len(moved):,} existing taxa to higher nodes.')

    logger.info('Reloading tree using new taxonomic positions.')
    tree.reload()

    logger.info('Identifying leaf nodes missing ranks')
    d_leaf_to_ranks_missing = tree.get_leaf_nodes_missing_rank()
    logger.info(f'Found {len(d_leaf_to_ranks_missing):,} genomes missing one or more ranks.')

    logger.info('Creating novel ranks')
    new_taxa_created = fill_in_missing_ranks(
        d_leaf_to_ranks_missing,
        tree.node_to_descs,
        min_bs,
        red_dict,
        metadata,
        set(tree.taxon_labels_to_node.keys())
    )
    logger.info(f'Created {len(new_taxa_created):,} novel ranks.')

    logger.info('Removing RED labels')
    tree.remove_red_labels()

    logger.info('Writing report to disk.')
    report = Report(moved, new_taxa_created)
    report.write(path_report_out)

    logger.info('Writing tree to disk.')
    tree.to_file(path_tree_out)

    logger.info('Done.')
    sys.exit(0)

    return



    shared = check_for_candidate_sharing(d_taxon_to_candidates)
    logger.info(f'There are {len(shared):,} taxa that share a candidate node.')

    moved = refine_nodes(d_taxon_to_candidates, tree.taxon_labels_to_node)
    moved = set(moved)
    logger.info(f'Moved {len(moved):,} taxa to higher nodes.')

    logger.info('Identifying leaf nodes missing ranks')
    d_leaf_to_ranks_missing = tree.get_leaf_nodes_missing_rank()

    logger.info('Filling in missing ranks')
    new_taxa_created = fill_in_missing_ranks(
        d_leaf_to_ranks_missing,
        tree.node_to_descs,
        min_bs,
        red_dict,
        metadata,
        set(tree.taxon_labels_to_node.keys())
    )
    logger.info('')




if __name__ == "__main__":
    app()
