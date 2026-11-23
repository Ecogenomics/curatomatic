"""
Limit a Rust report to the taxa that Rust and Python place differently.

Reads the Rust report (`<stem>.curatomatic.report.tsv`) and the output of
`scripts/compare_outputs.py` for the same dataset (`existing_taxa.tsv`, `novel_taxa.tsv`).

A taxon is placed the same when both versions put it on the same nodes: existing taxa in the
categories `unchanged` and `moved_same`, and novel taxa in `both_same_name` and
`both_different_name` (same node, the names differ only by Python's `NOVEL-` prefix, L4).

Writes to the output directory:
    report.placed_differently.tsv      the rows of the report placed differently, same columns.
                                       The report only lists taxa that Rust moved or created.
    report.placed_differently.all.tsv  every taxon placed differently, with the report columns plus
                                       the input, Python and Rust nodes and the likely cause. Taxa
                                       that Rust left alone have the action `not_moved` (Python
                                       moved them) or `python_only_novel` (only Python made them).

Only the standard library is used.

Usage:
    python scripts/report_placed_differently.py <report.tsv> <compare_outputs dir> <out dir>
"""
import argparse
import csv
from pathlib import Path

SAME = {'unchanged', 'moved_same', 'both_same_name', 'both_different_name'}
EXTRA = ['py_category', 'input_node', 'python_node', 'rust_node', 'cause']


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', restval='')
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('report', type=Path)
    parser.add_argument('compare_dir', type=Path)
    parser.add_argument('out_dir', type=Path)
    args = parser.parse_args()

    report = read_tsv(args.report)
    columns = list(report[0]) if report else []
    in_report = {r['taxon'] for r in report}
    existing = {r['taxon']: r for r in read_tsv(args.compare_dir / 'existing_taxa.tsv')}
    novel = read_tsv(args.compare_dir / 'novel_taxa.tsv')
    novel_rs = {r['rust']: r for r in novel if r['rust']}

    def placement(taxon):
        if taxon in existing:
            e = existing[taxon]
            return e['category'], e['input'], e['python'], e['rust'], e['cause']
        n = novel_rs[taxon]
        return n['category'], '', n['node'] if n['python'] else '', n['node'], n['cause']

    strict, full = [], []
    for row in report:
        values = placement(row['taxon'])
        if values[0] not in SAME:
            strict.append(row)
            full.append({**row, **dict(zip(EXTRA, values))})
    for taxon, e in existing.items():
        if e['category'] not in SAME and taxon not in in_report:
            values = (e['category'], e['input'], e['python'], e['rust'], e['cause'])
            full.append({'taxon': taxon, 'rank': e['rank'], 'action': 'not_moved', **dict(zip(EXTRA, values))})
    for n in novel:
        if n['category'] == 'python_only':
            cause = n['cause'] + (f" (rust has: {n['rust_existing']})" if n['rust_existing'] else '')
            full.append({'taxon': n['python'], 'rank': n['rank'], 'action': 'python_only_novel',
                         'n_genomes_after': n['n_genomes'], 'bootstrap_after': n['bootstrap'],
                         'red_after': n['red'], **dict(zip(EXTRA, (n['category'], '', n['node'], '', cause)))})

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir / 'report.placed_differently.tsv', strict, columns)
    write_tsv(args.out_dir / 'report.placed_differently.all.tsv', full, columns + EXTRA)
    print(f'{len(strict):,} report rows placed differently, {len(full):,} taxa placed differently in all')


if __name__ == '__main__':
    main()
