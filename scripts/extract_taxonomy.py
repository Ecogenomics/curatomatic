"""
Read the taxonomy of each genome off a decorated tree and write it to a TSV file.

For each leaf node whose label does not begin with "D-" (i.e. not a dummy node),
the taxon labels of every ancestor are collected and written as:

    <accession>\t<d__...; p__...; c__...; o__...; f__...; g__...; s__...>

Internal node labels may be in the form "SUPPORT:TAXA", "TAXA", or "SUPPORT".

Usage:
    python scripts/extract_taxonomy.py <tree.nwk> <output.tsv>
"""
import argparse
import sys
from pathlib import Path

import dendropy


def parse_taxa(label):
    """Return the list of taxa in an internal node label (empty if there are none)."""
    if not label:
        return []
    label = label.strip()
    if ':' in label:
        support, rest = label.split(':', 1)
        try:
            float(support)
            label = rest
        except ValueError:
            pass
    else:
        try:
            float(label)
            return []
        except ValueError:
            pass
    return [x.strip() for x in label.split(';') if x.strip()]


def extract_taxonomy(tree):
    """Yield (accession, taxonomy) for each non-dummy leaf in the tree."""
    cache = dict()

    def node_taxonomy(node):
        """The full taxonomy (root first) of an internal node, memoised."""
        if node is None:
            return []
        key = id(node)
        if key not in cache:
            cache[key] = node_taxonomy(node.parent_node) + parse_taxa(node.label)
        return cache[key]

    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100_000))
    for leaf in tree.leaf_node_iter():
        accession = leaf.taxon.label if leaf.taxon is not None else leaf.label
        if accession.startswith('D-'):
            continue
        yield accession, node_taxonomy(leaf.parent_node)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tree', type=Path, help='Path to the decorated Newick tree.')
    parser.add_argument('output', type=Path, help='Path to write the taxonomy TSV file.')
    args = parser.parse_args()

    tree = dendropy.Tree.get(
        path=str(args.tree),
        schema='newick',
        preserve_underscores=True,
    )

    n = 0
    with args.output.open('w') as fh:
        for accession, taxonomy in extract_taxonomy(tree):
            fh.write(f'{accession}\t{"; ".join(taxonomy)}\n')
            n += 1
    print(f'Wrote {n:,} genomes to {args.output}')


if __name__ == '__main__':
    main()
