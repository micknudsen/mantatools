from typing import Dict, Iterable, List

from svtoolbox.core import Variant
from svtoolbox.exceptions import InfoFieldNotFound


def parse_vcf(stream: Iterable) -> Dict[str, Variant]:
    """Read VCF file line by line and return a dictionary of with
    variant IDs as keys and Variant objects as values."""

    variants: Dict[str, Variant] = {}
    samples: List[str] = []

    for line in stream:

        # Skip header lines
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            continue

        columns = line.rstrip("\n").split("\t")

        try:
            format = columns[8]
            genoptypes = dict(zip(samples, columns[9:]))
        except IndexError:
            format = None
            genoptypes = None

        variant = Variant(
            chrom=columns[0],
            pos=columns[1],
            id=columns[2],
            ref=columns[3],
            alt=columns[4],
            qual=columns[5],
            filter=columns[6],
            info=columns[7],
            format=format,
            genotypes=genoptypes,
        )

        variants[variant.id] = variant

        # If the variant is a BND, and if we have already encountered
        # the mate variant, then link the two variants together.
        if variant.get_info("SVTYPE") == "BND":
            try:
                mate_id = str(variant.get_info("MATEID"))
                if mate_id in variants:
                    variant.mate = variants[mate_id]
                    variants[mate_id].mate = variant
            except InfoFieldNotFound:
                pass

    return variants
