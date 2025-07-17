import dendropy

from curatomatic.model.taxonomy import Taxonomy


class NodeLabel:

    def __init__(
            self,
            red: float,
            bootstrap: float | None,
            taxonomy: Taxonomy | None,
            gid: str | None
    ):
        self.red: float = red
        self.bootstrap: float | None = bootstrap
        self.taxonomy: Taxonomy | None = taxonomy
        self.gid: str | None = gid

    @classmethod
    def from_node(cls, node: dendropy.Node):
        if node.is_leaf():
            assert '|RED=' in node.taxon.label
            gid, red = node.taxon.label.strip().split('|RED=')
            red = float(red)
            assert red == 1.0
            return cls(red=red, bootstrap=None, taxonomy=None, gid=gid)
        else:
            assert '|RED=' in node.label
            prefix, suffix = node.label.strip().split('|RED=')
            red = float(suffix)
            bs = None
            taxonomy = None

            # Extract either the taxon or the bootstrap + taxon
            if ':' in prefix:
                bs, taxonomy = prefix.split(':')
                bs = float(bs)
            else:
                try:
                    bs = float(prefix)
                except ValueError:
                    taxonomy = prefix

            if taxonomy:
                taxonomy = Taxonomy.from_string(taxonomy)
            return cls(red=red, bootstrap=bs, taxonomy=taxonomy, gid=None)

    def to_string(self) -> str:
        tax_string = self.taxonomy.to_string() if self.taxonomy else ''
        if self.gid is not None:
            # This is a leaf node
            assert tax_string == ''
            return f'{self.gid}|RED={self.red:.3f}'
        else:
            # This is an internal node
            if self.bootstrap is None:
                assert tax_string != ''
                return f'{tax_string}|RED={self.red:.3f}'
            else:
                if tax_string == '':
                    return f'{self.bootstrap:.1f}|RED={self.red:.3f}'
                else:
                    return f'{self.bootstrap:.1f}:{tax_string}|RED={self.red:.3f}'
