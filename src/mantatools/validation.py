from typing import List

from pysam import AlignedSegment

from mantatools.core import Variant


def check_contig_support(variant: Variant, alignments: List[AlignedSegment]) -> bool:
    """Check if the variant is supported by alignment of contigs. This
    is very much work in progress and should not yet be used in production!"""

    # The contig does not even map to the refrence.
    if len(alignments) == 0:
        return False

    # The contig maps to the reference, and there is only the primary alignment.
    elif len(alignments) == 1:

        cigartuples = alignments[0].cigartuples or []

        match variant.get_info(key="SVTYPE"):

            case "DEL":
                return (2, abs(int(variant.get_info(key="SVLEN")))) in cigartuples

            case _:
                return False

    # The contig maps to the reference, and there is a least
    # one supplementary alignment.
    else:
        return True
