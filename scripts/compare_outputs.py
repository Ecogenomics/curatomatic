"""
Compare the output of the Python (v0.0.x) and Rust (v0.1.x) implementations of curatomatic.

Both versions are run on the same tree, so a node is identified by the set of genomes below it.
The key used is ``<lowest genome ID>/<number of genomes>``, which is unique for a tree without
unifurcations (two clades with the same lowest genome are nested, so they differ in size). Dummy
leaves (``D-``) are ignored, so the parent of a dummy pair has the same key as its genome: the
taxa it holds are the taxa of that genome's leaf.

Usage:
    python scripts/compare_outputs.py \\
        --input-tree   test_data/r226_ar53/gtdb_r226_ar53.decorated.scaled.red.tree \\
        --rs-input     test_data/r226_ar53/rs_input \\
        --py-out       test_data/r226_ar53/curatomatic_bs95.0 \\
        --rs-out       test_output/r226_ar53 \\
        --out          /tmp/compare/r226_ar53

``--input-tree`` is the original (old format) input tree, it gives the input taxa, RED and
bootstrap of each node. ``--rs-input`` is the converted input, it is only used as a sanity check.

Writes to ``--out``:
    summary.md          counts of each category, and the unexplained cases
    existing_taxa.tsv   each taxon of the input, and where each version put it
    novel_taxa.tsv      each novel taxon, matched between the versions by rank and node
    genome_diff.tsv     per rank, the (Python, Rust) taxon pairs of genomes that differ
    log_moves.tsv       the "Moving" lines of both logs, aligned per rank
"""
import argparse
import collections
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import dendropy

RANKS = 'dpcofgs'
RED_RE = re.compile(r'\|RED=[^|]*$', re.IGNORECASE)  # The converted Rust inputs use `red=`.
MOVE_RE = re.compile(r'Moving (\S+) from (.*) to (.*)$')
NOVEL_GID_RE = re.compile(r'^[a-z]__(NOVEL-)?G\d{9}$')


# --------------------------------------------------------------------------------------------------
# Taxa

def rank_of(taxon):
    return RANKS.index(taxon[0])


def is_placeholder(taxon):
    return len(taxon) == 3 and taxon.endswith('__')


def is_valid_name(taxon):
    """As ``Taxon::is_valid_name``: not a placeholder and not an unresolved species."""
    return not is_placeholder(taxon) and not taxon.startswith('s__{unresolved}')


def parse_label(label):
    """(bootstrap, [taxa]) of an internal node label, with or without a RED."""
    label = RED_RE.sub('', label or '').strip()
    bootstrap, taxa = None, label
    if ':' in label:
        bootstrap, taxa = label.split(':', 1)
    else:
        try:
            float(label)
            bootstrap, taxa = label, ''
        except ValueError:
            pass
    taxa = [x.strip() for x in taxa.split(';') if x.strip()]
    return (float(bootstrap) if bootstrap else None), [x for x in taxa if not is_placeholder(x)]


def parse_red(label):
    m = re.search(r'\|RED=([^|]*)$', label or '', re.IGNORECASE)
    return float(m.group(1)) if m else None


def gid_of(leaf):
    return RED_RE.sub('', leaf.taxon.label if leaf.taxon is not None else leaf.label)


# --------------------------------------------------------------------------------------------------
# Trees

@dataclass
class Node:
    key: str
    size: int
    parent: str | None
    taxa: list = field(default_factory=list)
    bootstrap: float | None = None
    red: float | None = None
    children: list = field(default_factory=list)


class Tree:
    """A tree reduced to its keyed nodes, with the taxa each node holds."""

    def __init__(self, path):
        sys.setrecursionlimit(1_000_000)
        tree = dendropy.Tree.get(path=str(path), schema='newick', preserve_underscores=True)
        self.nodes: dict[str, Node] = {}
        self.leaf_key: dict[str, str] = {}
        key_of = {}
        lo, size = {}, {}
        for n in tree.postorder_node_iter():
            if n.is_leaf():
                gid = gid_of(n)
                if gid.startswith('D-'):
                    continue
                lo[n], size[n] = gid, 1
            else:
                kids = [c for c in n.child_nodes() if c in lo]
                lo[n] = min(lo[c] for c in kids)
                size[n] = sum(size[c] for c in kids)
            key_of[n] = f'{lo[n]}/{size[n]}'

        for n in tree.preorder_node_iter():
            if n not in key_of:
                continue
            key = key_of[n]
            parent = n.parent_node
            while parent is not None and key_of[parent] == key:
                parent = parent.parent_node
            node = self.nodes.get(key)
            if node is None:
                node = Node(key, size[n], key_of[parent] if parent is not None else None)
                self.nodes[key] = node
                if node.parent is not None:
                    self.nodes[node.parent].children.append(key)
            if n.is_leaf():
                self.leaf_key[gid_of(n)] = key
            else:
                bootstrap, taxa = parse_label(n.label)
                node.taxa.extend(taxa)
                node.bootstrap = node.bootstrap if bootstrap is None else bootstrap
                red = parse_red(n.label)
                node.red = node.red if red is None else red

        # Euler tour, for ancestor tests.
        self.tin, self.tout = {}, {}
        root = next(k for k, v in self.nodes.items() if v.parent is None)
        t, stack = 0, [(root, False)]
        while stack:
            key, done = stack.pop()
            if done:
                self.tout[key] = t
                continue
            self.tin[key] = t
            t += 1
            stack.append((key, True))
            stack.extend((c, False) for c in reversed(self.nodes[key].children))
        self.root = root

    def add_leaf_taxa(self, taxonomy_path):
        """Put each row's taxa that are not on an ancestor onto its leaf (for trees without dummies)."""
        for gid, taxa in read_taxonomy(taxonomy_path).items():
            key = self.leaf_key[canonical(gid)]
            above = set()
            parent = self.nodes[key].parent
            while parent is not None:
                above.update(self.nodes[parent].taxa)
                parent = self.nodes[parent].parent
            self.nodes[key].taxa.extend(x for x in taxa if x not in above and not is_placeholder(x))

    def is_ancestor(self, a, b):
        """True if node ``a`` is a strict ancestor of ``b``."""
        return a != b and self.tin[a] <= self.tin[b] and self.tout[b] <= self.tout[a]

    def related(self, a, b):
        return a == b or self.is_ancestor(a, b) or self.is_ancestor(b, a)

    def placements(self):
        out = collections.defaultdict(list)
        for key in self.preorder():
            for taxon in self.nodes[key].taxa:
                out[taxon].append(key)
        return out

    def preorder(self):
        return sorted(self.nodes, key=self.tin.__getitem__)

    def genome_taxonomy(self):
        """The taxon nearest the leaf for each rank, as the Rust taxonomy file is written."""
        lineage = {self.root: self._apply({}, self.nodes[self.root].taxa)}
        for key in self.preorder():
            if key != self.root:
                lineage[key] = self._apply(lineage[self.nodes[key].parent], self.nodes[key].taxa)
        return {gid: lineage[key] for gid, key in self.leaf_key.items()}

    @staticmethod
    def _apply(parent, taxa):
        out = dict(parent)
        seen = set()
        for taxon in taxa:  # The first taxon of a rank on a node wins.
            r = taxon[0]
            if r not in seen:
                out[r] = taxon
                seen.add(r)
        return out


def canonical(gid):
    m = re.match(r'^(?:RS_|GB_)?GC[AF]_(\d{9})\.\d+$', gid)
    return f'G{m.group(1)}' if m else gid


def read_taxonomy(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.strip():
                gid, tax = line.rstrip('\n').split('\t')
                out[gid] = [x.strip() for x in tax.split(';') if x.strip()]
    return out


def read_report(path):
    with open(path, newline='') as fh:
        return {row['taxon']: row for row in csv.DictReader(fh, delimiter='\t')}


def read_moves(path):
    """The "Moving" lines of a log, per rank, in order."""
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            m = MOVE_RE.search(line.rstrip('\n'))
            if m:
                out[m.group(1)[0]].append((m.group(1), m.group(2), m.group(3)))
    return out


# --------------------------------------------------------------------------------------------------
# Comparison

def fmt_nodes(tree_in, keys):
    parts = []
    for k in keys:
        n = tree_in.nodes[k]
        red = '' if n.red is None else f'{n.red:.3f}'
        bs = '' if n.bootstrap is None else f'{n.bootstrap:g}'
        parts.append(f'{k}(red={red},bs={bs})')
    return ' '.join(parts)


def direction(tree, before, after):
    if before == after:
        return 'same'
    if len(before) != 1 or len(after) != 1:
        return 'multi'
    b, a = before[0], after[0]
    if tree.is_ancestor(a, b):
        return 'up'
    if tree.is_ancestor(b, a):
        return 'down'
    return 'across'


def nested_in(tree, keys, taxon):
    """The other named taxa of the same rank on an ancestor of any of ``keys``."""
    out = set()
    for key in keys:
        parent = tree.nodes[key].parent
        while parent is not None:
            out.update(x for x in tree.nodes[parent].taxa
                       if x[0] == taxon[0] and x != taxon and is_valid_name(x))
            parent = tree.nodes[parent].parent
    return sorted(out)


def valid_ranks_below(tree, exclude):
    """For each node, the ranks with a valid name at or below it (ignoring ``exclude``)."""
    out = {}
    for key in reversed(tree.preorder()):
        node = tree.nodes[key]
        ranks = {x[0] for x in node.taxa if is_valid_name(x) and x not in exclude}
        for c in node.children:
            ranks |= out[c]
        out[key] = ranks
    return out


def compare(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    t_in = Tree(args.input_tree)
    t_py = Tree(next(Path(args.py_out).glob('*.curatomatic.tree')))
    t_rs = Tree(next(Path(args.rs_out).glob('*.curatomatic.tree')))
    t_rs.add_leaf_taxa(next(Path(args.rs_out).glob('*.curatomatic.taxonomy.tsv')))
    t_rsin = Tree(next(Path(args.rs_input).glob('*.nwk')))
    t_rsin.add_leaf_taxa(next(Path(args.rs_input).glob('*.decorated_taxonomy.tsv')))
    report = read_report(next(Path(args.rs_out).glob('*.report.tsv')))

    sanity = []
    for name, t in (('python output', t_py), ('rust output', t_rs), ('rust input', t_rsin)):
        same = set(t.nodes) == set(t_in.nodes) and all(
            t.nodes[k].parent == t_in.nodes[k].parent for k in t_in.nodes)
        sanity.append(f'- Topology of the {name} matches the input: **{same}**')
    p_in, p_rsin = t_in.placements(), t_rsin.placements()
    sanity.append(f'- Taxa of the converted input match the original input: **{p_in == p_rsin}**')
    if p_in != p_rsin:
        diff = {k for k in set(p_in) | set(p_rsin) if p_in.get(k) != p_rsin.get(k)}
        sanity.append(f'  - {len(diff):,} differ, e.g. {sorted(diff)[:5]}')

    p_py, p_rs = t_py.placements(), t_rs.placements()
    in_taxa = set(p_in)
    novel_py = {x for x in p_py if x not in in_taxa}
    novel_rs = {x for x in p_rs if x not in in_taxa}

    # ---- Existing taxa.
    existing = []
    for taxon in sorted(in_taxa, key=lambda x: (rank_of(x), x)):
        if taxon[0] == 'd':
            continue
        ki, kp, kr = p_in[taxon], p_py.get(taxon, []), p_rs.get(taxon, [])
        si, sp, sr = set(ki), set(kp), set(kr)
        if not kp or not kr:
            cat = 'lost_in_python' if not kp else 'lost_in_rust'
            if not kp and not kr:
                cat = 'lost_in_both'
        elif sp == sr:
            cat = 'unchanged' if sp == si else 'moved_same'
        elif sr == si:
            cat = 'moved_python_only'
        elif sp == si:
            cat = 'moved_rust_only'
        else:
            cat = 'moved_differently'
        existing.append(dict(
            taxon=taxon, rank=taxon[0], category=cat,
            dir_python=direction(t_in, ki, kp) if kp else '',
            dir_rust=direction(t_in, ki, kr) if kr else '',
            input=fmt_nodes(t_in, ki), python=fmt_nodes(t_in, kp), rust=fmt_nodes(t_in, kr),
            rs_action=report.get(taxon, {}).get('action', ''),
            rs_flag=report.get(taxon, {}).get('flag', ''),
            nested_in=';'.join(nested_in(t_in, ki, taxon)),
            on_leaf_in_input='yes' if any(t_in.nodes[k].size == 1 for k in ki) else '',
        ))

    # ---- Novel taxa, matched by (rank, node).
    def novel_at(placements, taxa):
        return {(x[0], k): x for x in taxa for k in placements[x]}

    nv_py, nv_rs = novel_at(p_py, novel_py), novel_at(p_rs, novel_rs)
    below_rs = valid_ranks_below(t_rs, novel_rs)
    lineage_py, lineage_rs = t_py.genome_taxonomy(), t_rs.genome_taxonomy()
    leaves_under = collections.defaultdict(list)
    for gid, key in t_in.leaf_key.items():
        k = key
        while k is not None:
            leaves_under[k].append(gid)
            k = t_in.nodes[k].parent

    def overlapping(nv_other, rank, key):
        return sorted(f'{name}@{k}' for (r, k), name in nv_other.items()
                      if r == rank and k != key and t_in.related(k, key))

    def covered_by(lineage, rank, key, novel):
        """The existing (non-novel) taxa of ``rank`` that the genomes under ``key`` have."""
        return sorted({lineage[g].get(rank) for g in leaves_under[key]
                       if lineage[g].get(rank) and lineage[g].get(rank) not in novel})

    novel = []
    for rk in sorted(set(nv_py) | set(nv_rs), key=lambda x: (RANKS.index(x[0]), t_in.tin[x[1]])):
        rank, key = rk
        n_py, n_rs = nv_py.get(rk, ''), nv_rs.get(rk, '')
        if n_py and n_rs:
            cat = 'both_same_name' if n_py == n_rs else 'both_different_name'
        else:
            cat = 'python_only' if n_py else 'rust_only'
        lower = ''.join(r for r in RANKS[RANKS.index(rank) + 1:] if r in below_rs[key])
        node = t_in.nodes[key]
        novel.append(dict(
            rank=rank, node=key, n_genomes=node.size,
            red='' if node.red is None else f'{node.red:.3f}',
            bootstrap='' if node.bootstrap is None else f'{node.bootstrap:g}',
            python=n_py, rust=n_rs, category=cat,
            rs_flag=report.get(n_rs, {}).get('flag', '') if n_rs else '',
            valid_lower_ranks_below_in_rust=lower,
            other_python_novel=';'.join(overlapping(nv_py, rank, key)) if not n_py else '',
            other_rust_novel=';'.join(overlapping(nv_rs, rank, key)) if not n_rs else '',
            rust_existing=';'.join(covered_by(lineage_rs, rank, key, novel_rs)) if not n_rs else '',
            python_existing=';'.join(covered_by(lineage_py, rank, key, novel_py)) if not n_py else '',
        ))

    # ---- Attribution.
    diff_existing = [r for r in existing if r['category'] not in ('unchanged', 'moved_same')]
    diff_keys = {}  # taxon -> nodes involved, for cascade detection
    for r in diff_existing:
        t = r['taxon']
        diff_keys[t] = set(p_in[t]) | set(p_py.get(t, [])) | set(p_rs.get(t, []))
    for r in novel:
        if r['category'] != 'both_same_name':
            diff_keys[f'{r["rank"]}@{r["node"]}'] = {r['node']}

    def cascade_from(rank, keys, own):
        """Differences of a higher rank whose nodes are nested with ``keys``."""
        hits = []
        for name, other in diff_keys.items():
            if name == own or RANKS.index(name[0]) >= RANKS.index(rank):
                continue
            if any(t_in.related(a, b) for a in keys for b in other):
                hits.append(name)
        return sorted(hits)[:5]

    lost_rs = {r['taxon'] for r in existing if r['category'] == 'lost_in_rust'}

    def hidden_conflict(taxon):
        """Taxa of the same rank, below the Rust node of ``taxon``, that Rust never had."""
        return sorted(x for x in lost_rs if x[0] == taxon[0] and any(
            t_in.is_ancestor(a, b) for a in p_rs[taxon] for b in p_in[x]))

    for r in existing:
        cause = ''
        cat = r['category']
        if cat in ('unchanged', 'moved_same'):
            r['cause'] = ''
            r['possible_cascade_from'] = ''
            continue
        if r['rs_action'] == 'moved_down' or r['dir_rust'] == 'down':
            cause = 'L2 (downward refinement, Rust only)'
        elif cat in ('moved_python_only', 'moved_differently') and r['dir_python'] == 'down':
            cause = 'L2 (downward refinement, Python only)'
        elif cat == 'moved_python_only' and r['dir_python'] == 'up' and r['nested_in']:
            cause = 'L1 (nested taxon not moved up)'
        elif cat == 'lost_in_rust' and r['on_leaf_in_input']:
            cause = 'L6 (singleton not below the lowest rank on its path, not assigned to its leaf)'
        elif cat == 'moved_rust_only' and r['dir_rust'] == 'up' and hidden_conflict(r['taxon']):
            cause = ('L6 (a conflicting singleton below was not assigned to its leaf: '
                     + ';'.join(hidden_conflict(r['taxon'])) + ')')
        r['cause'] = cause
        r['possible_cascade_from'] = ';'.join(
            cascade_from(r['rank'], diff_keys[r['taxon']], r['taxon']))

    def missing_above_rs(rank, key):
        """The ranks above ``rank`` that a genome under ``key`` has no taxon for in Rust."""
        above = RANKS[1:RANKS.index(rank)]
        return ''.join(r for r in above if any(r not in lineage_rs[g] for g in leaves_under[key]))

    for r in novel:
        cause = ''
        cat = r['category']
        if cat == 'both_same_name':
            r['cause'] = r['possible_cascade_from'] = ''
            continue
        if cat == 'both_different_name':
            if NOVEL_GID_RE.match(r['python']) or r['rust'].endswith('__UNRESOLVABLE') \
                    or NOVEL_GID_RE.match(r['rust']):
                cause = 'L4 (novel name without a WGS prefix)'
        elif cat == 'python_only':
            if r['valid_lower_ranks_below_in_rust']:
                cause = 'L7 (novel taxon above a named taxon)'
            elif r['other_rust_novel']:
                cause = 'placed on another node by Rust'
            elif r['rust_existing']:
                cause = 'genomes hold an existing taxon in Rust'
            elif missing_above_rs(r['rank'], r['node']):
                cause = ('L7 (lower rank of a skipped novel taxon: Rust leaves '
                         + missing_above_rs(r['rank'], r['node']) + ' empty, so this rank is not filled)')
        elif cat == 'rust_only':
            if r['other_python_novel']:
                cause = 'placed on another node by Python'
            elif r['python_existing']:
                cause = 'genomes hold an existing taxon in Python'
        r['cause'] = cause
        r['possible_cascade_from'] = ';'.join(
            cascade_from(r['rank'], {r['node']}, f'{r["rank"]}@{r["node"]}'))

    # ---- Genome taxonomy differences.
    name_map = {r['python']: r['rust'] for r in novel if r['category'] == 'both_different_name'}
    genome_rows = collections.Counter()
    per_rank = {r: [0, 0] for r in RANKS}
    for gid in t_in.leaf_key:
        lp, lr = lineage_py[gid], lineage_rs[gid]
        for rank in RANKS:
            a, b = lp.get(rank, ''), lr.get(rank, '')
            if a != b:
                per_rank[rank][0] += 1
                mapped = name_map.get(a, a)
                if mapped != b:
                    per_rank[rank][1] += 1
                    genome_rows[(rank, a, b)] += 1

    # Sanity: the Rust taxonomy file is what the tree says.
    rs_file = read_taxonomy(next(Path(args.rs_out).glob('*.curatomatic.taxonomy.tsv')))
    n_bad = sum(1 for gid, taxa in rs_file.items()
                if [x for x in taxa if not is_placeholder(x)] !=
                [lineage_rs[canonical(gid)][r] for r in RANKS if r in lineage_rs[canonical(gid)]])
    sanity.append(f'- Rust taxonomy file matches the Rust tree: **{n_bad == 0}** ({n_bad:,} rows differ)')

    # ---- Logs.
    mv_py = read_moves(next(Path(args.py_out).glob('curatomatic.log')))
    mv_rs = read_moves(next(Path(args.rs_out).glob('curatomatic.log')))
    log_rows, first_div = [], {}
    for rank in RANKS:
        a, b = mv_py.get(rank, []), mv_rs.get(rank, [])
        for i in range(max(len(a), len(b))):
            x = a[i] if i < len(a) else None
            y = b[i] if i < len(b) else None
            if x != y and rank not in first_div:
                first_div[rank] = i
            log_rows.append(dict(rank=rank, index=i, same='yes' if x == y else '',
                                 python=' -> '.join(x) if x else '', rust=' -> '.join(y) if y else ''))

    # ---- Write.
    write_tsv(out / 'existing_taxa.tsv', existing)
    write_tsv(out / 'novel_taxa.tsv', novel)
    write_tsv(out / 'genome_diff.tsv', [
        dict(rank=r, python=a, rust=b, n_genomes=n)
        for (r, a, b), n in sorted(genome_rows.items(), key=lambda x: (RANKS.index(x[0][0]), -x[1]))])
    write_tsv(out / 'log_moves.tsv', log_rows)

    md = [f'# {args.name or out.name}: Python vs Rust', '', '## Sanity checks', '', *sanity, '']
    md += ['## Existing taxa', '', table(collections.Counter(r['category'] for r in existing)), '']
    md += ['### Causes of differences', '', table(collections.Counter(
        (r['category'], r['cause'] or 'UNEXPLAINED', 'cascade?' if r['possible_cascade_from'] else '')
        for r in existing if r['category'] not in ('unchanged', 'moved_same'))), '']
    md += ['## Novel taxa', '', table(collections.Counter((r['rank'], r['category']) for r in novel)), '']
    md += ['### Causes of differences', '', table(collections.Counter(
        (r['category'], r['cause'] or 'UNEXPLAINED', 'cascade?' if r['possible_cascade_from'] else '')
        for r in novel if r['category'] != 'both_same_name')), '']
    md += ['## Genomes with a different taxon', '',
           '| rank | genomes differ | after mapping novel names |', '| --- | --- | --- |',
           *[f'| {r} | {a:,} | {b:,} |' for r, (a, b) in per_rank.items()], '']
    md += ['## Refinement log ("Moving" lines)', '',
           '| rank | python | rust | first difference (index) |', '| --- | --- | --- | --- |',
           *[f'| {r} | {len(mv_py.get(r, []))} | {len(mv_rs.get(r, []))} | {first_div.get(r, "")} |'
             for r in 'pcofg'], '']
    (out / 'summary.md').write_text('\n'.join(md))
    print((out / 'summary.md').read_text())


def table(counter):
    rows = sorted(counter.items(), key=lambda x: (-x[1], str(x[0])))
    if not rows:
        return '_none_'
    width = len(rows[0][0]) if isinstance(rows[0][0], tuple) else 1
    lines = ['| ' + ' | '.join([''] * width) + ' | n |', '|' + ' --- |' * (width + 1)]
    for key, n in rows:
        key = key if isinstance(key, tuple) else (key,)
        lines.append('| ' + ' | '.join(map(str, key)) + f' | {n:,} |')
    return '\n'.join(lines)


def write_tsv(path, rows):
    with open(path, 'w', newline='') as fh:
        if not rows:
            return
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-tree', required=True, type=Path)
    parser.add_argument('--rs-input', required=True, type=Path)
    parser.add_argument('--py-out', required=True, type=Path)
    parser.add_argument('--rs-out', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--name', default=None)
    compare(parser.parse_args())


if __name__ == '__main__':
    main()
