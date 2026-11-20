"""
Convert the inputs of the Python curatomatic to the inputs of the Rust curatomatic.

For each release prefix (e.g. `r226_ar53`) in the input directory, the following are read:

    <prefix>.decorated.scaled.red.tree  - tree with "D-" dummy leaves
    <prefix>.decorated.dict             - JSON of the median RED of each rank
    <prefix>.metadata.tsv               - GTDB metadata with an ncbi_wgs_formatted column

and written to the output directory as:

    <prefix>.decorated.scaled.red.tree  - each dummy pair `('G|RED=1':0,'D-G|RED=1':0)'taxa|RED=1'`
                                          collapsed into the leaf `G`; internal labels unchanged
                                          apart from `; ` becoming `;` and `|RED=` becoming `|red=`
    <prefix>.taxonomy.tsv               - `<gid>\t<d__;p__;...;s__>` of each genome, from its path
                                          in the input tree (placeholders for missing ranks)
    <prefix>.median_red.tsv             - `rank\tmedian`, one row per rank from phylum to genus
    <prefix>.metadata.tsv               - only the accession and ncbi_wgs_formatted columns

Only the standard library is used. The tree is parsed iteratively, so it does not recurse.

Usage:
    python scripts/convert_python_inputs.py <input_dir> <output_dir>
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

RANKS = ('d', 'p', 'c', 'o', 'f', 'g', 's')
MEDIAN_RANKS = ('phylum', 'class', 'order', 'family', 'genus')

# A dummy pair and its labelled parent, the branch length after it is left in place.
DUMMY_RE = re.compile(
    r"\('(?P<gid>[^'|]+)\|RED=[^']*':[^,()]+,'D-(?P=gid)\|RED=[^']*':[^,()]+\)'(?P<label>[^']*)'"
)


def label_taxa(label):
    """The taxa of a label `[<bootstrap>][:][<taxa>]|RED=<red>`, highest rank first."""
    prefix = label.split('|', 1)[0]
    if ':' in prefix:
        prefix = prefix.split(':', 1)[1]
    else:
        try:
            float(prefix)
            return []
        except ValueError:
            pass
    return [t.strip() for t in prefix.split(';') if t.strip()]


def convert_label(label):
    """Normalise an internal label to the format of the Rust example trees."""
    if '|RED=' in label:
        prefix, red = label.split('|RED=', 1)
        prefix = ';'.join(t.strip() for t in prefix.split(';'))
        return f'{prefix}|red={red}'
    return label


def full_taxonomy(taxa):
    """Join taxa (highest rank first) into a seven rank string, with placeholders for gaps."""
    by_rank = {r: [] for r in RANKS}
    for taxon in taxa:
        by_rank[taxon[0]].append(taxon)
    out = []
    for rank in RANKS:
        out.extend(by_rank[rank] or [f'{rank}__'])
    return ';'.join(out)


def convert_tree(path_in, path_out):
    """Collapse the dummy pairs, and return the (gid, taxonomy) of each genome in tree order."""
    text = path_in.read_text().strip()

    leaf_taxa = {}

    def collapse(m):
        gid = m.group('gid')
        if gid in leaf_taxa:
            sys.exit(f'{path_in}: {gid} is given more than once')
        leaf_taxa[gid] = label_taxa(m.group('label'))
        return gid

    text = DUMMY_RE.sub(collapse, text)
    if "'D-" in text:
        sys.exit(f'{path_in}: a dummy leaf was not collapsed')

    # Scan the collapsed tree: track the taxa of each open clade, and assign them on close.
    out = []
    rows = []
    stack = [[]]  # the genomes below each open clade
    i, n = 0, len(text)
    ancestors_of = {}  # gid -> list of internal labels' taxa (lowest first)
    while i < n:
        ch = text[i]
        if ch == '(':
            stack.append([])
            out.append(ch)
            i += 1
        elif ch == ')':
            i += 1
            label = ''
            if i < n and text[i] == "'":
                j = text.index("'", i + 1)
                label = text[i + 1:j]
                i = j + 1
            taxa = label_taxa(label)
            gids = stack.pop()
            for gid in gids:
                ancestors_of[gid].append(taxa)
            stack[-1].extend(gids)
            out.append(')')
            if label:
                out.append(f"'{convert_label(label)}'")
        elif ch in ',:;':
            out.append(ch)
            i += 1
        else:
            # A leaf (genome ID) or a branch length.
            j = i
            while j < n and text[j] not in ',():;':
                j += 1
            token = text[i:j]
            if out and out[-1] == ':':
                out.append(token)
            else:
                if token not in leaf_taxa:
                    sys.exit(f'{path_in}: unexpected leaf {token!r}')
                ancestors_of[token] = []
                stack[-1].append(token)
                rows.append(token)
                out.append(token)
            i = j
    if len(stack) != 1:
        sys.exit(f'{path_in}: unbalanced parentheses')

    path_out.write_text(''.join(out) + '\n')

    result = []
    for gid in rows:
        taxa = []
        for anc in reversed(ancestors_of[gid]):
            taxa.extend(anc)
        taxa.extend(leaf_taxa[gid])
        result.append((gid, full_taxonomy(taxa)))
    return result


def convert_median(path_in, path_out):
    medians = json.loads(path_in.read_text())
    with path_out.open('w') as fh:
        fh.write('rank\tmedian\n')
        for rank in MEDIAN_RANKS:
            fh.write(f'{rank.capitalize()}\t{medians[rank]}\n')


def convert_metadata(path_in, path_out):
    csv.field_size_limit(sys.maxsize)
    with path_in.open(newline='') as fh_in, path_out.open('w') as fh_out:
        reader = csv.DictReader(fh_in, delimiter='\t')
        fh_out.write('accession\tncbi_wgs_formatted\n')
        for row in reader:
            fh_out.write(f"{row['accession']}\t{row['ncbi_wgs_formatted']}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_dir', type=Path)
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for tree in sorted(args.input_dir.glob('*.decorated.scaled.red.tree')):
        prefix = tree.name.removesuffix('.decorated.scaled.red.tree')
        print(f'Converting {prefix}')
        rows = convert_tree(tree, args.output_dir / tree.name)
        with (args.output_dir / f'{prefix}.taxonomy.tsv').open('w') as fh:
            for gid, tax in rows:
                fh.write(f'{gid}\t{tax}\n')
        convert_median(args.input_dir / f'{prefix}.decorated.dict', args.output_dir / f'{prefix}.median_red.tsv')
        convert_metadata(args.input_dir / f'{prefix}.metadata.tsv', args.output_dir / f'{prefix}.metadata.tsv')
        print(f'  {len(rows):,} genomes')


if __name__ == '__main__':
    main()
