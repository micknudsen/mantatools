from typing import List

from pysam import AlignedSegment

from mantatools.core import Variant


def check_contig_support(variant: Variant, alignments: List[AlignedSegment]) -> bool:
    """Check if the variant is supported by alignment of contigs. This
    is very much work in progress and should not yet be used in production!"""
    return len(alignments) > 0
