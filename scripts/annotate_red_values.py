"""
This script will transfer RED values from the node_rd.tsv output from Phylorank
on to an existing tree.

This is generally not needed for running curatomatic, although this may be helpful for debugging.
"""
from pathlib import Path

import dendropy
from tqdm import tqdm

from curatomatic.model.node_rd import NodeRd
from curatomatic.model.tree import Tree
from curatomatic.util.tree import mrca_faster


def mrca_to_red(node_rd, tree):
    changes = dict()
    d_taxon_to_node = {x.taxon.label: x for x in tree.leaf_node_iter()}

    d_depth_to_nodes = Tree.get_depth_to_nodes(tree)
    d_node_desc_nodes = Tree.get_node_desc_nodes(d_depth_to_nodes)

    for gid_a, gid_b, red in tqdm(node_rd.data):

        if gid_a == gid_b:
            # Just annotate the Taxon with this
            changes[d_taxon_to_node[gid_a]] = red

        else:
            # Annotate internal nodes
            mrca = mrca_faster(gid_a, gid_b, d_taxon_to_node, d_node_desc_nodes, tree)
            if mrca in changes:
                if changes[mrca] != red:
                    print('???')
            changes[mrca] = red

    return changes


def annotate(d_mrca_to_red):
    next_id = 0
    for mrca, red in tqdm(d_mrca_to_red.items()):
        if mrca.is_leaf():
            mrca.taxon.label += f'_xx_RED={red}'
        else:
            if mrca.label is None:
                mrca.label = 's__MISSING ' + f'uid{next_id}_xx_RED={red}'
                next_id += 1
            else:
                mrca.label += f'_xx_RED={red}'
    return

def fix_quoted(path):
    # Newer versions of dendropy don't quote leaf nodes correctly
    with path.open() as f:
        content = f.read()
    content = content.replace('_xx_', '|')
    with path.open('w') as f:
        f.write(content)

def main():
    domain = 'ar53'
    root_dir = Path(f'/Users/aaron/projects/curatomatic/data/r226/{domain}')
    path_tree_in = root_dir / Path(f'gtdb_r226_{domain}.decorated.scaled.tree')
    path_node_rd = root_dir / Path(f'gtdb_r226_{domain}.decorated.node_rd.tsv')
    path_tree_out = root_dir / Path(f'gtdb_r226_{domain}.decorated.scaled.red.tree')

    print('Reading node rd file')
    node_rd = NodeRd.from_file(path_node_rd)

    print('Reading tree')
    tree = dendropy.Tree.get(path=str(path_tree_in), schema="newick", preserve_underscores=True)
    tree.is_rooted = True

    print('Finding MRCA nodes')
    d_mrca_to_red = mrca_to_red(node_rd, tree)
    annotate(d_mrca_to_red)

    tree.write(path=str(path_tree_out), schema="newick", suppress_rooting=True, unquoted_underscores=False)
    fix_quoted(path_tree_out)

    return


if __name__ == '__main__':
    main()
