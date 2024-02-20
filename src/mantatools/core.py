from dataclasses import dataclass
from typing import Dict, List, Union

from mantatools.exceptions import GenotypeFieldNotFound, InfoFieldNotFound


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
        self.genotype_dicts: List[Dict[str, str]] = []

        for entry in self.info.split(";"):
            if "=" in entry:
                key, value = entry.split("=", 1)
                self.info_dict[key] = value
            else:
                self.info_dict[entry] = True

        for genotype in self.genotypes:
            genotype_dict = {}
            for i, key in enumerate(self.format.split(":")):
                genotype_dict[key] = genotype.split(":")[i]
            self.genotype_dicts.append(genotype_dict)

    def get_info(self, key: str) -> Union[str, bool]:
        try:
            return self.info_dict[key]
        except KeyError:
            raise InfoFieldNotFound()

    def get_genotype(self, key: str, sample: int = 0) -> str:
        try:
            return self.genotype_dicts[sample][key]
        except KeyError:
            raise GenotypeFieldNotFound()
