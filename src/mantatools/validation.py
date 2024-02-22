from typing import List

from pysam import AlignedSegment

from mantatools.core import Variant


def check_contig_support(variant: Variant, alignments: List[AlignedSegment]) -> bool:
    return True
