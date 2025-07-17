from pathlib import Path


class NodeRd:

    def __init__(self, data):
        self.data = data

    @classmethod
    def from_file(cls, path: Path):
        out = list()
        with path.open() as f:
            for line in f.readlines():
                gid, red = line.strip().split('\t')
                red = float(red)
                if '|' in gid:
                    taxon_a, taxon_b = gid.split('|')
                else:
                    taxon_a = gid
                    taxon_b = gid
                out.append((taxon_a, taxon_b, red))
        return cls(data=out)
