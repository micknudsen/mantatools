from dataclasses import dataclass
from typing import Dict, List, Union

from mantatools.exceptions import InfoFieldNotFound


@dataclass
class Variant:
    chrom: str
    pos: str
    id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: str
    format: str
    genotypes: List[str]

    def __post_init__(self) -> None:
        self.info_dict: Dict[str, Union[str, bool]] = {}
        for entry in self.info.split(";"):
            if "=" in entry:
                key, value = entry.split("=", 1)
                self.info_dict[key] = value
            else:
                self.info_dict[entry] = True

    def get_info(self, key: str) -> Union[str, bool]:
        try:
            return self.info_dict[key]
        except KeyError:
            raise InfoFieldNotFound()
