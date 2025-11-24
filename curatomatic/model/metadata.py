from collections import defaultdict
from pathlib import Path

from curatomatic.model.taxonomy import Taxonomy
from curatomatic.util.genome import gid_to_canonical


class MetadataRow:
    __slots__ = ('accession', 'taxonomy', 'wgs')

    def __init__(self, accession: str, taxonomy: Taxonomy | None, wgs: str | None):
        self.accession: str = accession
        self.taxonomy: Taxonomy | None = taxonomy
        self.wgs: str | None = wgs


class Metadata:

    def __init__(self, rows: dict[str, MetadataRow]):
        self.rows: dict[str, MetadataRow] = rows

    @classmethod
    def from_file(cls, path: Path):
        out = dict()
        with path.open() as f:
            header = f.readline().strip().split('\t')
            header = {k: i for i, k in enumerate(header)}
            for line in f.readlines():
                cols = line.strip().split('\t')
                gid = gid_to_canonical(cols[header['accession']])

                taxonomy, wgs_formatted = None, None
                if len(cols) > 1:
                    if 'ncbi_wgs_formatted' in header:
                        wgs_formatted = cols[header['ncbi_wgs_formatted']]
                        wgs_formatted = None if wgs_formatted == 'none' else wgs_formatted

                    if 'gtdb_taxonomy' in header:
                        if cols[header['gtdb_taxonomy']] != 'none':
                            taxonomy = Taxonomy.from_string(cols[header['gtdb_taxonomy']])

                out[gid] = MetadataRow(
                    accession=gid,
                    taxonomy=taxonomy,
                    wgs=wgs_formatted
                )
        return cls(rows=out)

    def subset_to_ids(self, ids: set[str]):
        self.rows = {k: v for k, v in self.rows.items() if k in ids}

    def get_taxonomy_from_gids(self, gids: set[str]) -> dict[str, Taxonomy]:
        out = dict()
        for gid in gids:
            if gid in self.rows:
                cur_row = self.rows[gid]
                if cur_row.taxonomy is not None:
                    out[gid] = cur_row.taxonomy
        return out

    def get_taxon_to_gids(self, gids: set[str]) -> dict[str, set[str]]:
        out = defaultdict(set)
        for gid in gids:
            if gid in self.rows:
                cur_row = self.rows[gid]
                if cur_row.taxonomy is not None:
                    for taxon in cur_row.taxonomy.iter_ranks():
                        out[taxon].add(gid)
        return dict(out)