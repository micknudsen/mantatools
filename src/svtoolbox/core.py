from dataclasses import dataclass
from typing import Dict, List, Optional, Union

from svtoolbox.exceptions import (
    FieldNotFound,
    GenotypeFieldNotFound,
    InfoFieldNotFound,
    MissingMate,
)


@dataclass
class Position:
    chrom: str
    pos: int

    def __str__(self) -> str:
        return f"{self.chrom}:{self.pos}"


@dataclass
class Interval:
    chrom: str
    left: int
    right: int

    def overlaps(self, other: "Interval") -> bool:
        """VCF intervals are 1-based and closed"""
        return (
            self.chrom == other.chrom
            and self.right >= other.left
            and self.left <= other.right
        )

    def __str__(self) -> str:
        return f"{self.chrom}:{self.left}-{self.right}"


@dataclass
class BedPE:
    chrom_1: str
    start_1: int
    end_1: int
    chrom_2: str
    start_2: int
    end_2: int
    name: str
    score: Optional[str]
    strand_1: Optional[str]
    strand_2: Optional[str]
    fields: Optional[Dict[str, str]] = None

    @classmethod
    def from_intervals(
        cls,
        left: Interval,
        right: Interval,
        name: str,
        score: Optional[str],
        strand_1: Optional[str],
        strand_2: Optional[str],
        fields: Optional[Dict[str, str]] = None,
    ) -> "BedPE":
        """Convenience method for creating a BedPE object from two intervals.
        Note that the BEDPE format is 0-based and half-open, whereas the VCF
        format is 1-based and closed."""
        return cls(
            chrom_1=left.chrom,
            start_1=left.left - 1,
            end_1=left.right,
            chrom_2=right.chrom,
            start_2=right.left - 1,
            end_2=right.right,
            name=name,
            score=score,
            strand_1=strand_1,
            strand_2=strand_2,
            fields=fields,
        )

    def __str__(self) -> str:
        columns = [
            self.chrom_1,
            str(self.start_1),
            str(self.end_1),
            self.chrom_2,
            str(self.start_2),
            str(self.end_2),
            self.name,
            self.score if self.score is not None else ".",
            self.strand_1 if self.strand_1 is not None else ".",
            self.strand_2 if self.strand_2 is not None else ".",
        ]
        if self.fields is not None:
            columns.append(
                ";".join([f"{key}={value}" for key, value in self.fields.items()])
            )
        return "\t".join(columns)


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
    format: Optional[str]
    genotypes: Optional[Dict[str, str]]

    # The mate variant of a BND variant
    mate: Optional["Variant"] = None

    def __post_init__(self) -> None:
        """Gather INFO and FORMAT values in dictionaries. This is convenient
        for looking values up later. There is one dictionary for INFO fields
        and one FORMAT dictionary for each sample in the VCF file."""

        self.info_dict: Dict[str, Union[str, bool]] = {}
        self.format_dicts: Dict[str, Dict[str, str]] = {}

        for entry in self.info.split(";"):
            # This is a key-value entry
            if "=" in entry:
                key, value = entry.split("=", 1)
                self.info_dict[key] = value
            # This is a flag entry
            else:
                self.info_dict[entry] = True

        if self.format is not None and self.genotypes is not None:
            for sample, genotypes in self.genotypes.items():
                self.format_dicts[sample] = {
                    key: value
                    for key, value in zip(self.format.split(":"), genotypes.split(":"))
                }

    def __str__(self) -> str:
        rows: List[str] = [
            self.chrom,
            self.pos,
            self.id,
            self.ref,
            self.alt,
            self.qual,
            self.filter,
            ";".join(
                [
                    f"{key}={value}" if isinstance(value, str) else key
                    for key, value in self.info_dict.items()
                ]
            ),
        ]

        if self.format is not None and self.genotypes is not None:
            rows.append(self.format)
            rows.extend(self.genotypes.values())

        return "\t".join(rows)

    def get_info(self, key: str) -> Union[str, bool]:
        """Return the value of the INFO field with the given key. If the key
        is a flag, return True. If the key is not found, raise an exception."""
        try:
            return self.info_dict[key]
        except KeyError:
            raise InfoFieldNotFound(key)

    def set_info(self, key: str, value: Union[str, bool]) -> None:
        """Set the value of the INFO field with the given key.."""
        self.info_dict[key] = value

    def get_genotype(
        self,
        sample: str,
        key: str,
    ) -> str:
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
            # This is the style used by Delly
            try:
                chrom, pos = self.get_info("CHR2"), self.get_info("POS2")
                return Position(str(chrom), int(pos))
            except InfoFieldNotFound:
                pass
            # This is the style used by Manta
            if self.mate is None:
                raise MissingMate(self.id)
            return self.mate.start
        return Position(self.chrom, int(self.get_info("END")))

    @property
    def ci_start(self) -> Interval:
        """Return confidence interval for the start position. If the CIPOS
        info field is not found, return an interval with the start position
        as both left and right."""
        try:
            left, right = str(self.get_info("CIPOS")).split(",")
            return Interval(
                chrom=self.chrom,
                left=self.start.pos + int(left),
                right=self.start.pos + int(right),
            )
        except InfoFieldNotFound:
            return Interval(
                chrom=self.chrom,
                left=self.start.pos,
                right=self.start.pos,
            )

    @property
    def ci_end(self) -> Interval:
        """Return confidence interval for the end position. If the CIEND
        info field is not found, return an interval with the end position
        as both left and right. For BND variants, the confidence interval
        of the end poistion is the same as the confidence interval of the
        start position of the mate."""
        if self.get_info("SVTYPE") == "BND" and self.mate is not None:
            return self.mate.ci_start
        try:
            left, right = str(self.get_info("CIEND")).split(",")
            return Interval(
                chrom=self.end.chrom,
                left=self.end.pos + int(left),
                right=self.end.pos + int(right),
            )
        except InfoFieldNotFound:
            return Interval(
                chrom=self.end.chrom,
                left=self.end.pos,
                right=self.end.pos,
            )

    def to_bedpe(self, include_fields: Optional[List[str]] = None) -> BedPE:
        """Create a BEDPE representation of the variant."""

        fields: Dict[str, str] = {}
        if include_fields is not None:
            for name in include_fields:
                match name:
                    case "REF":
                        fields[name] = self.ref
                    case "ALT":
                        fields[name] = self.alt
                    case "QUAL":
                        fields[name] = self.qual
                    case "FILTER":
                        fields[name] = self.filter
                    case _:
                        raise FieldNotFound(name)

        return BedPE.from_intervals(
            left=self.ci_start,
            right=self.ci_end,
            name=self.id,
            score=self.qual,
            strand_1=None,
            strand_2=None,
            fields=fields if fields else None,
        )
