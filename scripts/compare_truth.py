"""
Score the pre-curation of curatomatic (Rust, and optionally Python) against a curated taxonomy.

The curated taxonomy (the "truth") is what the curators published after they fixed the
pre-curation by hand. The input taxonomy is scored too, as a baseline ("no curatomatic at all").
Only the genomes on the tree are scored.

For each genome and rank (d to g), the taxon nearest the leaf is compared with the truth:

    correct   same name as the truth
    wrong     a different name
    empty     no name (a placeholder): the truth never has one, so this is always a miss

Novel taxa without a WGS name are ``<rank>__NOVEL-<gid>`` in Python and ``<rank>__<gid>`` in Rust.
Python's are renamed to Rust's form before scoring.

Names alone can't credit a novel taxon that a curator renamed, so each rank is also scored by its
grouping, ignoring names:

    pair precision   of the pairs of genomes sharing a taxon in the output, the fraction that
                     share one in the truth
    pair recall      of the pairs of genomes sharing a taxon in the truth, the fraction that
                     share one in the output (empty genomes share nothing)
    exact clades     output taxa whose genome set is exactly a truth taxon's (on the tree)

Species are scored separately: curatomatic never changes them, so they measure the input.

Usage:
    python scripts/compare_truth.py \\
        --name     r226_ar53 \\
        --truth    test_data/r226_ar53/ar53_taxonomy_r226_curated.tsv \\
        --baseline test_data/r226_ar53/rs_input/gtdb_r226_ar53.decorated_taxonomy.tsv \\
        --rs-out   <rust output dir> \\
        --py-out   test_data/r226_ar53/curatomatic_bs95.0 \\
        --out      <out dir>/r226_ar53

``--baseline`` is a taxonomy file (literal path or one taxon per rank): the last taxon of each
rank is the one nearest the leaf. ``--py-out`` is optional. No dependencies beyond the standard library.

Writes to ``--out``:
    scores.tsv     per version and rank: counts, pair precision and recall, exact clades
    genomes.tsv    per genome and rank where any version differs from the truth
    taxa.tsv       per version, rank and output taxon: its best matching truth taxon
    changes.tsv    per version and rank, genomes fixed or broken relative to the baseline,
                   grouped by (baseline, output, truth) taxon
    py_vs_rs.tsv   per rank, genomes where exactly one of Python and Rust is correct, grouped by
                   (Python, Rust, truth) taxon
"""
import argparse
import collections
import csv
import re
from pathlib import Path

RANKS = 'dpcofg'
NOVEL_RE = re.compile(r'^([a-z])__NOVEL-(G\d{9})$')


def canonical(gid):
    m = re.match(r'^(?:RS_|GB_)?GC[AF]_(\d{9})\.\d+$', gid)
    return f'G{m.group(1)}' if m else gid


def is_placeholder(taxon):
    return len(taxon) == 3 and taxon.endswith('__')


def read_rows(path):
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.strip():
                gid, tax = line.split('\t')[:2]
                yield canonical(gid), [x.strip() for x in tax.split(';') if x.strip()]


def nearest_leaf(taxa):
    """{rank: taxon} with the last named taxon of each rank (placeholders dropped)."""
    out = {}
    for taxon in taxa:
        if not is_placeholder(taxon):
            out[taxon[0]] = taxon
    return out


def read_taxonomy(path):
    return {gid: nearest_leaf(taxa) for gid, taxa in read_rows(path)}


def node_taxa(label):
    """The named taxa of an internal node label: ``[bootstrap:]taxa[|RED=x]``, any case of RED."""
    label = re.sub(r'\|red=[^|]*$', '', label, flags=re.IGNORECASE).strip()
    if ':' in label:
        label = label.split(':', 1)[1]
    elif re.fullmatch(r'[\d.]+', label):
        return []
    return [x.strip() for x in label.split(';') if x.strip() and not is_placeholder(x.strip())]


def newick_genome_taxonomy(path):
    """{gid: {rank: taxon}} of a Newick tree, as ``compare_outputs.Tree.genome_taxonomy``.

    For each rank, the taxon nearest the leaf; on a node with several taxa of a rank, the first.
    Dummy leaves (``D-``) are skipped: the parent of a dummy pair holds its genome's leaf taxa.
    """
    text = Path(path).read_text().strip()
    out = {}
    leaves = []  # (gid, index of its parent clade)
    clades = []  # label of each clade, by index, in the order they open
    open_ids = []  # indices of the clades open at the current position
    parents = {}
    i, n = 0, len(text)
    next_id = 0

    def read_label(i):
        if i < n and text[i] == "'":
            j, buf = i + 1, []
            while True:
                k = text.index("'", j)
                buf.append(text[j:k])
                if k + 1 < n and text[k + 1] == "'":
                    buf.append("'")
                    j = k + 2
                else:
                    i = k + 1
                    break
            label = ''.join(buf)
        else:
            j = i
            while j < n and text[j] not in ':,();':
                j += 1
            label, i = text[i:j], j
        if i < n and text[i] == ':':
            i += 1
            while i < n and text[i] not in ',();':
                i += 1
        return label, i

    while i < n:
        c = text[i]
        if c == '(':
            parents[next_id] = open_ids[-1] if open_ids else None
            open_ids.append(next_id)
            clades.append('')
            next_id += 1
            i += 1
        elif c == ',':
            i += 1
        elif c == ')':
            label, i = read_label(i + 1)
            clades[open_ids.pop()] = label
        elif c == ';':
            break
        else:
            label, i = read_label(i)
            leaves.append((label, open_ids[-1]))

    lineage = {}
    for cid in range(len(clades)):  # parents are always opened before their children
        parent = parents[cid]
        lin = dict(lineage[parent]) if parent is not None else {}
        seen = set()
        for taxon in node_taxa(clades[cid]):
            if taxon[0] not in seen:
                lin[taxon[0]] = taxon
                seen.add(taxon[0])
        lineage[cid] = lin
    for gid, parent in leaves:
        gid = re.sub(r'\|red=[^|]*$', '', gid, flags=re.IGNORECASE)
        if not gid.startswith('D-'):
            out[canonical(gid)] = lineage[parent]
    return out


def read_python(py_out):
    taxa = newick_genome_taxonomy(next(Path(py_out).glob('*.curatomatic.tree')))
    return {gid: {r: NOVEL_RE.sub(r'\1__\2', x) for r, x in t.items()} for gid, t in taxa.items()}


# --------------------------------------------------------------------------------------------------
# Scores

def pairs(n):
    return n * (n - 1) // 2


def score_rank(gids, out, truth, rank):
    counts = collections.Counter()
    cells = collections.Counter()
    out_groups = collections.defaultdict(set)
    truth_groups = collections.defaultdict(set)
    for gid in gids:
        t = truth[gid][rank]
        truth_groups[t].add(gid)
        o = out[gid].get(rank)
        if o is None:
            counts['empty'] += 1
            continue
        counts['correct' if o == t else 'wrong'] += 1
        out_groups[o].add(gid)
        cells[(o, t)] += 1

    tp = sum(pairs(n) for n in cells.values())
    out_pairs = sum(pairs(len(g)) for g in out_groups.values())
    truth_pairs = sum(pairs(len(g)) for g in truth_groups.values())
    truth_sets = {frozenset(g) for g in truth_groups.values()}
    exact = sum(1 for g in out_groups.values() if frozenset(g) in truth_sets)
    return {
        'genomes': len(gids),
        'correct': counts['correct'],
        'wrong': counts['wrong'],
        'empty': counts['empty'],
        'pct_correct': pct(counts['correct'], len(gids)),
        'pair_precision': pct(tp, out_pairs),
        'pair_recall': pct(tp, truth_pairs),
        'taxa_out': len(out_groups),
        'taxa_truth': len(truth_groups),
        'exact_clades': exact,
        'pct_exact_clades': pct(exact, len(out_groups)),
    }, out_groups, truth_groups


def pct(a, b):
    return f'{100 * a / b:.3f}' if b else ''


def best_matches(version, rank, out_groups, truth_groups, input_taxa):
    truth_of = {g: t for t, gs in truth_groups.items() for g in gs}
    rows = []
    for taxon, gs in sorted(out_groups.items()):
        overlap = collections.Counter(truth_of[g] for g in gs)
        best, n = overlap.most_common(1)[0]
        union = len(gs | truth_groups[best])
        rows.append({
            'version': version, 'rank': rank, 'taxon': taxon,
            'novel': taxon not in input_taxa,
            'n_genomes': len(gs), 'best_truth': best, 'n_truth': len(truth_groups[best]),
            'overlap': n, 'jaccard': f'{n / union:.4f}',
            'exact_clade': gs == truth_groups[best], 'same_name': taxon == best,
            'n_truth_taxa': len(overlap),
        })
    return rows


def grouped(rows, limit_examples=3):
    """Group (key tuple, gid) rows, most common first, with example genomes."""
    groups = collections.defaultdict(list)
    for k, gid in rows:
        groups[k].append(gid)
    return sorted(((k, len(v), ','.join(sorted(v)[:limit_examples])) for k, v in groups.items()),
                  key=lambda x: (-x[1], x[0]))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--name', required=True)
    p.add_argument('--truth', required=True, type=Path)
    p.add_argument('--baseline', required=True, type=Path)
    p.add_argument('--rs-out', required=True, type=Path)
    p.add_argument('--py-out', type=Path)
    p.add_argument('--out', required=True, type=Path)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rs_path = next(args.rs_out.glob('*.curatomatic.taxonomy.tsv'))
    versions = {'baseline': read_taxonomy(args.baseline), 'rust': read_taxonomy(rs_path)}
    if args.py_out:
        versions['python'] = read_python(args.py_out)

    gids = sorted(versions['rust'])
    for name, tax in versions.items():
        if set(tax) != set(gids):
            raise SystemExit(f'{name}: genomes differ from the Rust output')
    all_truth = dict(read_rows(args.truth))
    missing = [g for g in gids if g not in all_truth]
    if missing:
        raise SystemExit(f'{len(missing):,} tree genomes not in the truth, e.g. {missing[:3]}')
    truth = {}
    for gid in gids:
        t = nearest_leaf(all_truth[gid])
        if set(t) != set(RANKS + 's') or t['s'].startswith('s__{unresolved}'):
            raise SystemExit(f'{gid}: the truth is incomplete: {all_truth[gid]}')
        truth[gid] = t
    input_taxa = {x for t in versions['baseline'].values() for x in t.values()}

    scores, taxa = [], []
    for version, tax in versions.items():
        for rank in RANKS + 's':
            s, out_groups, truth_groups = score_rank(gids, tax, truth, rank)
            scores.append({'dataset': args.name, 'version': version, 'rank': rank, **s})
            if rank in 'pcofg' and version != 'baseline':
                taxa.extend(best_matches(version, rank, out_groups, truth_groups, input_taxa))
    write_tsv(args.out / 'scores.tsv', scores)
    write_tsv(args.out / 'taxa.tsv', taxa)

    names = list(versions)
    genome_rows = []
    for gid in gids:
        for rank in RANKS:
            vals = {v: versions[v][gid].get(rank, '') for v in names}
            if any(x != truth[gid][rank] for x in vals.values()):
                genome_rows.append({'gid': gid, 'rank': rank, 'truth': truth[gid][rank], **vals})
    write_tsv(args.out / 'genomes.tsv', genome_rows, ['gid', 'rank', 'truth', *names])

    change_rows = []
    for version in names[1:]:
        for rank in RANKS:
            rows = []
            for gid in gids:
                b, o, t = (versions['baseline'][gid].get(rank, ''), versions[version][gid].get(rank, ''),
                           truth[gid][rank])
                if (b == t) != (o == t):
                    rows.append(((('fixed' if o == t else 'broken'), b, o, t), gid))
            for (kind, b, o, t), n, ex in grouped(rows):
                change_rows.append({'version': version, 'rank': rank, 'change': kind, 'baseline': b,
                                    'output': o, 'truth': t, 'n_genomes': n, 'examples': ex})
    write_tsv(args.out / 'changes.tsv', change_rows,
              ['version', 'rank', 'change', 'baseline', 'output', 'truth', 'n_genomes', 'examples'])

    if 'python' in versions:
        rows = []
        for rank in RANKS:
            diffs = []
            for gid in gids:
                py, rs, t = versions['python'][gid].get(rank, ''), versions['rust'][gid].get(rank, ''), \
                    truth[gid][rank]
                if (py == t) != (rs == t):
                    diffs.append((('python' if py == t else 'rust', py, rs, t), gid))
            for (winner, py, rs, t), n, ex in grouped(diffs):
                rows.append({'rank': rank, 'correct': winner, 'python': py, 'rust': rs, 'truth': t,
                             'n_genomes': n, 'examples': ex})
        write_tsv(args.out / 'py_vs_rs.tsv', rows,
                  ['rank', 'correct', 'python', 'rust', 'truth', 'n_genomes', 'examples'])

    print(f'{args.name}: {len(gids):,} genomes scored, written to {args.out}')


if __name__ == '__main__':
    main()
