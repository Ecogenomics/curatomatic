from collections import defaultdict, deque

import dendropy

from curatomatic.model.metadata import Metadata
from curatomatic.model.node_label import NodeLabel
from curatomatic.model.rank import RANKS
from curatomatic.model.red import RED
from curatomatic.model.taxonomy import Taxonomy
from curatomatic.util.taxon import generate_unique_rank


def parse_node_label(node: dendropy.Node):
    red = None
    bs = None
    taxon = None

    # Leaf node, return nothing
    if node.is_leaf():
        red = 1.0
        bs = None
        taxon = None

    # Dummy node for leaf nodes, extract taxon and return RED
    else:
        label = node.label

        # Extract the RED that will always be present
        prefix, suffix = label.split('|RED=')
        red = float(suffix)

        # Extract either the taxon or the bootstrap + taxon
        if ':' in prefix:
            bs, taxon = prefix.split(':')
            bs = float(bs)
        else:
            try:
                bs = float(prefix)
            except ValueError:
                taxon = prefix

    return red, bs, taxon


def get_taxa_for_node(node, d_node_labels):
    out = defaultdict(set)
    red, bs, taxa = d_node_labels[node]
    if taxa:
        for taxon in taxa.split('; '):
            out[taxon[0]].add(taxon)
    return out


def check_for_candidate_sharing(d_taxon_to_candidates):
    out = dict()
    ranks = ['p', 'c', 'o', 'f', 'g', 's']
    for rank in ranks:
        seen = set()
        shared = set()

        for taxon, node_lst in d_taxon_to_candidates.items():
            if taxon[0] != rank:
                continue

            node_set = set(node_lst)
            common = seen.intersection(node_set)
            if len(common) > 0:
                shared.update(common)
            seen.update(node_lst)

        for taxon, node_lst in d_taxon_to_candidates.items():
            if taxon[0] != rank:
                continue

            node_set = set(node_lst)
            common = shared.intersection(node_set)
            if len(common) > 0:
                out[taxon] = common
    return out


def refine_nodes(
        d_taxon_to_candidates: dict[str, list[dendropy.Node]],
        d_taxon_to_node: dict[str, dendropy.Node]
):
    n_moved = 0
    moved = list()
    # This needs to be done from the highest ranks first to avoid any conflicts
    for cur_rank in 'pcofgs':
        for taxon, lst_candidates in d_taxon_to_candidates.items():
            if taxon[0] == cur_rank:
                new_node = lst_candidates[-1]
                old_node = d_taxon_to_node[taxon]
                transfer_taxon(old_node, new_node, taxon)
                moved.append(taxon)
                n_moved += 1
    return moved


def transfer_taxon(old_node: dendropy.Node, new_node: dendropy.Node, taxon: str):
    old_info = NodeLabel.from_node(old_node)
    new_info = NodeLabel.from_node(new_node)

    if new_info.taxonomy is None:
        new_info.taxonomy = Taxonomy()

    old_info.taxonomy.remove_taxon(taxon)
    new_info.taxonomy.add_taxon(taxon)

    old_node.label = old_info.to_string()
    new_node.label = new_info.to_string()
    return


def get_taxa_set_for_desc_nodes(d_node_to_descs, cur_node):
    out = defaultdict(list)
    for desc_node in d_node_to_descs[cur_node]:
        desc_tax = NodeLabel.from_node(desc_node).taxonomy
        _, _, tax = parse_node_label(desc_node)
        if desc_tax:
            for cur_taxon in desc_tax.iter_ranks():
                out[cur_taxon[0]].append(cur_taxon)
    return out


def fill_in_missing_ranks(d_leaf_to_ranks_missing, d_node_to_descs, min_bs, red_dict: RED, metadata: Metadata,
                          known_taxa):
    d_rank_to_leaf_candidates = defaultdict(dict)

    for leaf_node, missing_ranks in d_leaf_to_ranks_missing.items():

        # Try find the highest node placement for each rank descending
        for cur_rank in sorted(missing_ranks, key=lambda x: RANKS.index(x)):

            # Go up the tree to find a suitable node that wont conflict
            candidates = list()
            cur_node = leaf_node.parent_node
            while cur_node is not None:
                cur_desc_taxa = get_taxa_set_for_desc_nodes(d_node_to_descs, cur_node)
                if len(cur_desc_taxa[cur_rank]) > 0:
                    break
                else:
                    # Check RED and bootstral support
                    cur_node_label = NodeLabel.from_node(cur_node)
                    cur_bs = cur_node_label.bootstrap
                    cur_red = cur_node_label.red

                    bs_supported = cur_bs is None or cur_bs >= min_bs

                    parent_red = red_dict.data[RANKS[RANKS.index(cur_rank) - 1]]
                    midpoint_red = parent_red + ((red_dict.data[cur_rank] - parent_red) / 2)
                    red_supported = cur_red >= midpoint_red

                    if bs_supported and red_supported:
                        candidates.append(cur_node)
                cur_node = cur_node.parent_node

            # Store the candidate nodes
            d_rank_to_leaf_candidates[cur_rank][leaf_node] = candidates

    # Process higher to lower ranks
    new_taxa_created = set()
    for rank in sorted(d_rank_to_leaf_candidates.keys(), key=lambda x: RANKS.index(x)):

        # Check to see if any other taxa share the same candidate node
        d_candidate_to_leaves = defaultdict(set)
        for leaf_node, candidates in d_rank_to_leaf_candidates[rank].items():
            d_candidate_to_leaves[candidates[-1]].add(leaf_node)

        # Create novel groups for that highest candidate node that covers desc taxa
        # Frist, we need to generate a name
        for candidate, leaf_node_set in d_candidate_to_leaves.items():
            leaf_node_gids = {NodeLabel.from_node(x).gid for x in leaf_node_set}
            cur_known_taxa = known_taxa.union(new_taxa_created)
            new_taxon = generate_unique_rank(rank, leaf_node_gids, metadata, cur_known_taxa)
            new_taxa_created.add(new_taxon)

            candidate_info = NodeLabel.from_node(candidate)

            # Format the label
            if not candidate_info.taxonomy:
                candidate_info.taxonomy = Taxonomy()
            candidate_info.taxonomy.add_taxon(new_taxon)

            # Add the node label
            candidate.label = candidate_info.to_string()

    return new_taxa_created


def create_novel_ranks(tree, red_dict, d_leaf_to_ranks_missing):
    # Group the missing ranks by the rank they are missing
    d_rank_to_leaf_nodes_missing = defaultdict(set)
    for leaf_node, missing_ranks in d_leaf_to_ranks_missing.items():
        for rank in missing_ranks:
            d_rank_to_leaf_nodes_missing[rank].add(leaf_node)

    # Find the highest nodes that we can possibly use for each of the cases
    # starting from p to s
    for cur_rank, missing_leaf_nodes in sorted(d_rank_to_leaf_nodes_missing.items(), key=lambda x: RANKS.index(x[0])):
        print()

        # Find the highest node that we can use for each of the missing ranks
        d_taxon_to_candidates = defaultdict(list)

    return


def mrca_faster(
        gid_a: str,
        gid_b: str,
        d_gid_to_node: dict[str, dendropy.Node],
        d_node_desc_nodes: dict[dendropy.Node, set[dendropy.Node]],
        tree: dendropy.Tree
) -> dendropy.Node:
    node_a = d_gid_to_node[gid_a]
    node_b = d_gid_to_node[gid_b]

    queue = deque([tree.seed_node])

    last_node = None
    while len(queue) > 0:
        cur_node = queue.popleft()
        last_node = cur_node

        for child_node in cur_node.child_nodes():
            child_descs = d_node_desc_nodes[child_node]

            # If they both exist down this path, keep going
            if node_a in child_descs and node_b in child_descs:
                queue.append(child_node)
    return last_node


def get_tax_strings_for_leaf_nodes(node: dendropy.Node, d_gid_to_taxonomy) -> set[str]:
    out = set()



    return out