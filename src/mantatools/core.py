from dataclasses import dataclass
from typing import List


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
