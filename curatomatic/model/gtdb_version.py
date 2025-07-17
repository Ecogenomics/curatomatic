from enum import Enum
from pathlib import Path


class GtdbVersion(str, Enum):
    R207 = 'R207'
    R214 = 'R214'
    R220 = 'R220'
    R226_precuration = 'R226-precuration'

    def get_metadata_path(self) -> Path:
        return Path(f'curatomatic/data/gtdb/{self.value}_metadata.tsv.gz')
