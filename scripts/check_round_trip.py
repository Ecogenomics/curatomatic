"""
Check that a tree with dummy nodes is unchanged after a round trip through curatomatic.

The round trip is: remove the dummy nodes of an old style tree (`scripts/convert_python_inputs.py`),
load the tree and its taxonomy file with `curatomatic --debug`, and read the tree it writes with
the dummy nodes put back (`<stem>.debug.zombie.tree`). The two trees must be the same: the same
topology, child order, branch lengths, bootstraps, RED values and taxa on every node, and the same
dummy pair above every genome.

Labels are compared by their values, not their text: `; ` and `;` both separate taxa, `RED=` and
`red=` are the same key, and `1.0` and `1.000` are the same RED. Branch lengths are compared as
text, since curatomatic never changes them.

Optionally, two taxonomy files can be compared too (e.g. the input `--taxonomy` file and
`<stem>.debug.taxonomy.tsv`): they must have the same rows, in any order.

Only the standard library is used. The trees are parsed iteratively, so it does not recurse.

Usage:
    python scripts/check_round_trip.py <original.tree> <debug.zombie.tree> [<taxonomy.tsv> <debug.taxonomy.tsv>]
"""
import argparse
import sys
from pathlib import Path


def tokenize(text):
    """Yield the tokens of a Newick string: `(`, `,`, `)`, `;`, `:<length>` and labels."""
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in '(),;':
            yield c
            i += 1
        elif c == ':':
            j = i + 1
            while j < n and text[j] not in '(),;':
                j += 1
            yield ':' + text[i + 1:j].strip()
            i = j
        elif c == "'":
            j = i + 1
            label = []
            while True:
                k = text.index("'", j)
                label.append(text[j:k])
                if k + 1 < n and text[k + 1] == "'":
                    label.append("'")
                    j = k + 2
                else:
                    break
            yield 'L' + ''.join(label)
            i = k + 1
        elif c.isspace():
            i += 1
        else:
            j = i
            while j < n and text[j] not in "(),;:":
                j += 1
            yield 'L' + text[i:j].strip().replace('_', ' ')
            i = j


def parse(path):
    """Parse a tree into nodes `[label, length, children]`, returning the root."""
    root = [None, None, []]
    stack = [root]
    last = None
    for tok in tokenize(Path(path).read_text()):
        if tok == '(':
            node = [None, None, []]
            stack[-1][2].append(node)
            stack.append(node)
            last = None
        elif tok == ',':
            last = None
        elif tok == ')':
            last = stack.pop()
        elif tok == ';':
            break
        elif tok.startswith(':'):
            (last if last is not None else stack[-1][2][-1])[1] = tok[1:]
        else:
            if last is None:
                # A leaf.
                node = [tok[1:], None, []]
                stack[-1][2].append(node)
                last = node
            else:
                last[0] = tok[1:]
    # The outermost brackets are the root.
    if len(root[2]) != 1:
        sys.exit(f'{path}: expected one root')
    return root[2][0]


def label_value(label):
    """A label as (bootstrap, taxa, red), each None if it is not given."""
    if label is None:
        return (None, (), None)
    red = None
    for key in ('|RED=', '|red=', ':RED=', ':red='):
        if key in label:
            label, red = label.split(key, 1)
            red = float(red)
            if key.startswith(':'):
                label += ':'
            break
    bootstrap = None
    if ':' in label:
        bootstrap, label = label.split(':', 1)
    else:
        try:
            float(label)
            bootstrap, label = label, ''
        except ValueError:
            pass
    if bootstrap is not None and bootstrap != '':
        bootstrap = float(bootstrap)
    elif bootstrap == '':
        bootstrap = None
    taxa = tuple(t.strip() for t in label.split(';') if t.strip())
    return (bootstrap, taxa, red)


def walk(root):
    """Yield (path, node) for every node in pre-order, without recursion."""
    stack = [('0', root)]
    while stack:
        path, node = stack.pop()
        yield path, node
        for i, child in reversed(list(enumerate(node[2]))):
            stack.append((f'{path}.{i}', child))


def compare_trees(path_a, path_b, max_report=10):
    a, b = parse(path_a), parse(path_b)
    n_nodes = n_leaves = n_dummy = 0
    problems = []
    stack = [('root', a, b)]
    while stack:
        where, x, y = stack.pop()
        n_nodes += 1
        if len(x[2]) != len(y[2]):
            problems.append(f'{where}: {len(x[2])} children vs {len(y[2])}')
            continue
        if x[1] != y[1]:
            problems.append(f'{where}: branch length {x[1]!r} vs {y[1]!r}')
        if not x[2]:
            n_leaves += 1
            n_dummy += x[0].startswith('D-')
            if label_value(x[0]) != label_value(y[0]):
                problems.append(f'{where}: leaf {x[0]!r} vs {y[0]!r}')
        else:
            vx, vy = label_value(x[0]), label_value(y[0])
            if vx != vy:
                problems.append(f'{where}: label {x[0]!r} vs {y[0]!r}')
            for i, (cx, cy) in enumerate(zip(x[2], y[2])):
                stack.append((f'{where}.{i}', cx, cy))
    print(f'Trees: {n_nodes:,} nodes, {n_leaves:,} leaves ({n_dummy:,} dummy), {len(problems):,} differences')
    for p in problems[:max_report]:
        print(f'  {p}')
    return not problems


def read_rows(path):
    rows = {}
    for line in Path(path).read_text().splitlines():
        if line.strip():
            gid, tax = line.split('\t')
            rows[gid] = ';'.join(t.strip() for t in tax.split(';'))
    return rows


def compare_taxonomy(path_a, path_b, max_report=10):
    a, b = read_rows(path_a), read_rows(path_b)
    only_a = sorted(set(a) - set(b))
    only_b = sorted(set(b) - set(a))
    diff = sorted(g for g in set(a) & set(b) if a[g] != b[g])
    print(f'Taxonomy: {len(a):,} vs {len(b):,} rows, {len(only_a):,} only in the first, '
          f'{len(only_b):,} only in the second, {len(diff):,} different')
    for g in diff[:max_report]:
        print(f'  {g}\n    {a[g]}\n    {b[g]}')
    return not (only_a or only_b or diff)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('original_tree', type=Path)
    parser.add_argument('zombie_tree', type=Path)
    parser.add_argument('taxonomy', type=Path, nargs='?')
    parser.add_argument('debug_taxonomy', type=Path, nargs='?')
    args = parser.parse_args()

    ok = compare_trees(args.original_tree, args.zombie_tree)
    if args.taxonomy and args.debug_taxonomy:
        ok = compare_taxonomy(args.taxonomy, args.debug_taxonomy) and ok
    print('OK' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
