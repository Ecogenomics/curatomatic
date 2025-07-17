from collections import defaultdict


class Taxonomy:

    __slots__ = ('d', 'p', 'c', 'o', 'f', 'g', 's')

    def __init__(
            self,
            d: list[str] | None = None,
            p: list[str] | None = None,
            c: list[str] | None = None,
            o: list[str] | None = None,
            f: list[str] | None = None,
            g: list[str] | None = None,
            s: list[str] | None = None
    ):
        self.d = d if d else list()
        self.p = p if p else list()
        self.c = c if c else list()
        self.o = o if o else list()
        self.f = f if f else list()
        self.g = g if g else list()
        self.s = s if s else list()

    def __repr__(self):
        return self.to_string()

    @classmethod
    def from_string(cls, tax_string: str):
        out = defaultdict(list)
        taxa = tax_string.split(';')
        for taxon in taxa:
            taxon = taxon.strip()
            out[taxon[0]].append(taxon)

        if any([len(x) > 1 for x in out.values()]):
            print(f'Warning: Shared name: {taxa}')
            # raise Exception('Not implemented.')

        return cls(
            d=out['d'],
            p=out['p'],
            c=out['c'],
            o=out['o'],
            f=out['f'],
            g=out['g'],
            s=out['s'],
        )

    def iter_ranks(self):
        for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
            val = getattr(self, rank)
            for taxon in val:
                yield taxon

    def remove_taxon(self, taxon: str):
        taxon_rank = taxon[0]
        cur_val = getattr(self, taxon_rank)
        cur_val.remove(taxon)

    def add_taxon(self, taxon: str):
        taxon_rank = taxon[0]
        cur_val = getattr(self, taxon_rank)
        cur_val.append(taxon)

    def to_string(self):
        out = list()
        for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
            cur_val = getattr(self, rank)
            out.extend(cur_val)
        return '; '.join(out)

    def all_assigned(self) -> bool:
        for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
            cur_val = getattr(self, rank)
            if len(cur_val) == 0:
                return False
        return True

    def missing_ranks(self) -> set[str]:
        out = set()
        for rank in ('d', 'p', 'c', 'o', 'f', 'g', 's'):
            cur_val = getattr(self, rank)
            if len(cur_val) == 0:
                out.add(rank)
        return out

    def get_rank(self, rank: str) -> list[str]:
        return getattr(self, rank)