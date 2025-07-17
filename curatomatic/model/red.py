import json
from pathlib import Path

from curatomatic.model.rank import RANKS


class RED:

    def __init__(self, data):
        self.data: dict[str, float] = data

    @classmethod
    def from_file(cls, path: Path):
        with path.open() as f:
            val = json.load(f)
        out = {k[0]: v for k, v in val.items()}
        out['d'] = 0.0
        return cls(data=out)

    def to_str(self):
        out = list()
        for k, v in sorted(self.data.items(), key=lambda x: RANKS.index(x[0])):
            out.append(f'{k}={v:.3f}')
        return ', '.join(out)