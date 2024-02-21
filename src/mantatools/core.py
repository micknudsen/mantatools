from dataclasses import dataclass
from typing import Dict, List, Optional, Union

from mantatools.exceptions import GenotypeFieldNotFound, InfoFieldNotFound, MissingMate


@dataclass
class Position:
    chrom: str
    pos: int


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

    # The mate variant of a BND variant
    mate: Optional["Variant"] = None

    def __post_init__(self) -> None:
        """Gather INFO and FORMAT values in dictionaries. This is convenient
        for looking values up later. There is one dictionary for INFO fields
        and one FORMAT dictionary for each sample in the VCF file."""

        self.info_dict: Dict[str, Union[str, bool]] = {}
        self.format_dicts: List[Dict[str, str]] = []

        for entry in self.info.split(";"):
            # This is a key-value entry
            if "=" in entry:
                key, value = entry.split("=", 1)
                self.info_dict[key] = value
            # This is a flag entry
            else:
                self.info_dict[entry] = True

        for genotype in self.genotypes:
            self.format_dicts.append(
                {
                    key: value
                    for key, value in zip(self.format.split(":"), genotype.split(":"))
                }
            )

    def get_info(self, key: str) -> Union[str, bool]:
        """Return the value of the INFO field with the given key. If the key
        is a flag, return True. If the key is not found, raise an exception."""
        try:
            return self.info_dict[key]
        except KeyError:
            raise InfoFieldNotFound(key)

    def get_genotype(self, key: str, sample: int = 0) -> str:
        """Return the value of the FORMAT field with the given key for the
        given sample. If the key is not found, raise an exception."""
        try:
            return self.format_dicts[sample][key]
        except KeyError:
            raise GenotypeFieldNotFound(key)

    @property
    def start(self) -> Position:
        """Return the start postion of the variant."""
        return Position(self.chrom, int(self.pos))

    @property
    def end(self) -> Position:
        """Return the end position of the variant. For BND variants, this
        is the start position of the mate. For all other variants, the
        end position is specified in the END info field."""
        if self.get_info("SVTYPE") == "BND":
            if self.mate is None:
                raise MissingMate(self.id)
            return self.mate.start
        return Position(self.chrom, int(self.get_info("END")))

    @property
    def ci_start(self) -> Interval:
        raise NotImplementedError()
