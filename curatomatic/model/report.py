from pathlib import Path

from curatomatic.model.rank import RANKS


class Report:

    def __init__(self, moved, new_taxa_created):
        self.moved = moved
        self.new_taxa_created = new_taxa_created

    def write(self, path: Path):
        with path.open('w') as f:
            for taxon in sorted(self.moved.union(self.new_taxa_created), key=lambda x: RANKS.index(x[0])):
                if taxon in self.moved:
                    f.write(f'{taxon}\tMOVED\n')
                else:
                    f.write(f'{taxon}\tCREATED\n')
