"""
All tests in this file are based on real variants called by Manta. FASTQ files
were downloaded (https://doi.org/10.1101/2020.12.11.422022), and the following
steps were performed:

1. Sequencing adapters were removed using cutadapt
2. Reads were mapped to the hg38 reference genome using bwa mem
3. Duplicate reads were marked using samtools markdup
4. Variants were called using Manta

https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG001.novaseq.pcr-free.30x.R1.fastq.gz
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG001.novaseq.pcr-free.30x.R2.fastq.gz

The idea is to use test-driven development using a manually curated set of
variants to create (almost) perfect algorithms for assesing structural variants.
"""

import unittest

from pysam import AlignedSegment, AlignmentHeader

from svtoolbox.core import Variant
from svtoolbox.validation import check_contig_support


class TestContigSupport(unittest.TestCase):
    def setUp(self) -> None:
        self.header = AlignmentHeader.from_dict(
            {
                "SQ": [
                    {"SN": "chr1", "LN": 248956422},
                    {"SN": "chr2", "LN": 242193529},
                    {"SN": "chr3", "LN": 198295559},
                    {"SN": "chr4", "LN": 190214555},
                    {"SN": "chr5", "LN": 181538259},
                    {"SN": "chr6", "LN": 170805979},
                    {"SN": "chr7", "LN": 159345973},
                    {"SN": "chr8", "LN": 145138636},
                    {"SN": "chr9", "LN": 138394717},
                    {"SN": "chr10", "LN": 133797422},
                    {"SN": "chr11", "LN": 135086622},
                    {"SN": "chr12", "LN": 133275309},
                    {"SN": "chr13", "LN": 114364328},
                    {"SN": "chr14", "LN": 107043718},
                    {"SN": "chr15", "LN": 101991189},
                    {"SN": "chr16", "LN": 90338345},
                    {"SN": "chr17", "LN": 83257441},
                    {"SN": "chr18", "LN": 80373285},
                    {"SN": "chr19", "LN": 58617616},
                    {"SN": "chr20", "LN": 64444167},
                    {"SN": "chr21", "LN": 46709983},
                    {"SN": "chr22", "LN": 50818468},
                    {"SN": "chrX", "LN": 156040895},
                    {"SN": "chrY", "LN": 57227415},
                    {"SN": "chrM", "LN": 16569},
                ]
            }
        )

    def test_check_contig_support_true_1(self) -> None:
        """This is an example of a clear-cut case, where the contig supports the
        variant. The contig remaps as two soft-clipped reads, one at the start
        and one at the end of the contig, with the positions of the soft-clipped
        bases exactly matching the breakpoints of the variant."""

        variant = Variant(
            chrom="chr18",
            pos="57279506",
            id="MantaDEL:51284:0:1:0:0:0",
            ref="A",
            alt="<DEL>",
            qual="999",
            filter="PASS",
            info="END=57281486;SVTYPE=DEL;SVLEN=-1980;CONTIG=AATAGGGATTGACATAGCTCTGGCAGAAAAGTCCTGGGGAAGTTTCTCATTGGCTTAGCTGGAGGTAAACTCCCAGCTCTGCACTAATCACTGGGCCAGAGGCATGGCTGTTTGGACAAAAAGGAACAGGTTCTGCTCATTATATCTAGGAAGATCCTACTTTGTTTCAGGCACTGGGCTGACATCTATTCTCTCTTATCTTTTCGACATTGAGAAGTAGACATCATGATCTGTCCTTTATAGATGAAGAGACAAAGGCCAACACTCCTCTCCCAACACTGGACAT",
            format="GT:FT:GQ:PL:PR:SR",
            genotypes=["1/1:PASS:77:999,80,0:0,26:2,19"],
        )

        alignments = [
            AlignedSegment.from_dict(
                {
                    "name": "MantaDEL:51284:0:1:0:0:0",
                    "flag": "0",
                    "ref_name": "chr18",
                    "ref_pos": "57281487",
                    "map_quality": "60",
                    "cigar": "139S147M",
                    "next_ref_name": "*",
                    "next_ref_pos": "0",
                    "length": "0",
                    "seq": "AATAGGGATTGACATAGCTCTGGCAGAAAAGTCCTGGGGAAGTTTCTCATTGGCTTAGCTGGAGGTAAACTCCCAGCTCTGCACTAATCACTGGGCCAGAGGCATGGCTGTTTGGACAAAAAGGAACAGGTTCTGCTCATTATATCTAGGAAGATCCTACTTTGTTTCAGGCACTGGGCTGACATCTATTCTCTCTTATCTTTTCGACATTGAGAAGTAGACATCATGATCTGTCCTTTATAGATGAAGAGACAAAGGCCAACACTCCTCTCCCAACACTGGACAT",
                    "qual": "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                    "tags": [
                        "NM:i:0",
                        "ms:i:294",
                        "AS:i:294",
                        "nn:i:0",
                        "tp:A:P",
                        "cm:i:20",
                        "s1:i:137",
                        "s2:i:0",
                        "de:f:0",
                        "SA:Z:chr18,57279368,+,139M147S,60,0;",
                        "rl:i:0",
                    ],
                },
                header=self.header,
            ),
            AlignedSegment.from_dict(
                {
                    "name": "MantaDEL:51284:0:1:0:0:0",
                    "flag": "2048",
                    "ref_name": "chr18",
                    "ref_pos": "57279368",
                    "map_quality": "60",
                    "cigar": "139M147S",
                    "next_ref_name": "*",
                    "next_ref_pos": "0",
                    "length": "0",
                    "seq": "AATAGGGATTGACATAGCTCTGGCAGAAAAGTCCTGGGGAAGTTTCTCATTGGCTTAGCTGGAGGTAAACTCCCAGCTCTGCACTAATCACTGGGCCAGAGGCATGGCTGTTTGGACAAAAAGGAACAGGTTCTGCTCATTATATCTAGGAAGATCCTACTTTGTTTCAGGCACTGGGCTGACATCTATTCTCTCTTATCTTTTCGACATTGAGAAGTAGACATCATGATCTGTCCTTTATAGATGAAGAGACAAAGGCCAACACTCCTCTCCCAACACTGGACAT",
                    "qual": "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                    "tags": [
                        "NM:i:0",
                        "ms:i:278",
                        "AS:i:278",
                        "nn:i:0",
                        "tp:A:P",
                        "cm:i:14",
                        "s1:i:132",
                        "s2:i:0",
                        "de:f:0",
                        "SA:Z:chr18,57281487,+,139S147M,60,0;",
                        "rl:i:0",
                    ],
                },
                header=self.header,
            ),
        ]

        self.assertTrue(
            check_contig_support(
                variant=variant,
                alignments=alignments,
            )
        )

    def test_check_contig_support_true_2(self) -> None:
        """This is an example of a clear-cut case, where the contig supports the
        variant. The contig remaps as a single which spans the entire variant."""

        variant = Variant(
            chrom="chr4",
            pos="1092989",
            id="MantaDEL:31576:0:0:0:0:0",
            ref="CCAGGGTCTGCGTGCTCAGTGCTGACGCAGCCTGTGGTAGGGCAGAGGCTG",
            alt="C",
            qual="999",
            filter="PASS",
            info="END=1093039;SVTYPE=DEL;SVLEN=-50;CIGAR=1M50D;CONTIG=ACTCACGCGGACTCTCGCCAAGAGGCCAGGAGAGGCGGCTGCCTGGTCCGGAGCACACTTCTCACTCTTCGGTTCAATACCAGTTCTCCTCCATGGAGTGGCCTGTGCCTGCATTCGTCCACATGAGCTCCCGACTACGCCAGGACGGTGGATGGAACGGATCTTAGAGGATTACTGGGAAGAGGAAGACGTTTAATTGTTACCAACTAGACTAGGAAGGTGGACCAGGGCAGGCAGTGGCTGAGGTGGTTTGTTGACGCTGTCCCAGGGCAGGTCCTGAGGCCTGAGC",
            format="GT:FT:GQ:PL:PR:SR",
            genotypes=["1/1:PASS:108:999,111,0:0,0:0,41"],
        )

        alignments = [
            AlignedSegment.from_dict(
                {
                    "name": "MantaDEL:31576:0:0:0:0:0",
                    "flag": "0",
                    "ref_name": "chr4",
                    "ref_pos": "1092844",
                    "map_quality": "60",
                    "cigar": "146M50D143M",
                    "next_ref_name": "*",
                    "next_ref_pos": "0",
                    "length": "0",
                    "seq": "ACTCACGCGGACTCTCGCCAAGAGGCCAGGAGAGGCGGCTGCCTGGTCCGGAGCACACTTCTCACTCTTCGGTTCAATACCAGTTCTCCTCCATGGAGTGGCCTGTGCCTGCATTCGTCCACATGAGCTCCCGACTACGCCAGGACGGTGGATGGAACGGATCTTAGAGGATTACTGGGAAGAGGAAGACGTTTAATTGTTACCAACTAGACTAGGAAGGTGGACCAGGGCAGGCAGTGGCTGAGGTGGTTTGTTGACGCTGTCCCAGGGCAGGTCCTGAGGCCTGAGC",
                    "qual": "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                    "tags": [
                        "NM:i:53",
                        "ms:i:534",
                        "AS:i:474",
                        "nn:i:0",
                        "tp:A:P",
                        "cm:i:31",
                        "s1:i:220",
                        "s2:i:0",
                        "de:f:0.0138",
                        "rl:i:0",
                    ],
                },
                header=self.header,
            ),
        ]

        self.assertTrue(
            check_contig_support(
                variant=variant,
                alignments=alignments,
            )
        )

    def test_check_contig_support_false_1(self) -> None:
        """This is an example of an IMPRECISE variant which by definitions does not
        have a contig seqeunce. Thus is variant is obviously not supported by one."""

        variant = Variant(
            chrom="chr13",
            pos="16000650",
            id="MantaDEL:10:24582:24583:0:0:0",
            ref="A",
            alt="<DEL>",
            qual="39",
            filter="PASS",
            info="END=16006192;SVTYPE=DEL;SVLEN=-5542;IMPRECISE;CIPOS=-230,231;CIEND=-147,148",
            format="GT:FT:GQ:PL:PR",
            genotypes=["0/1:PASS:39:89,0,47:4,5"],
        )

        alignments = []

        self.assertFalse(
            check_contig_support(
                variant=variant,
                alignments=alignments,
            )
        )

    def test_check_contig_support_false_2(self) -> None:
        """This is an example of a variant not supported by the contig. The primary alignment
        does not span the variant, and there is no secondary alignment."""

        variant = Variant(
            chrom="chr1",
            pos="34102573",
            id="MantaDEL:2026:0:0:0:2:0",
            ref="GATCCAGCCATATCCCTGTAATCCGGCCATATCCCTGTGATCCAACCATATCCCTGTAATCCGGCCATATCCCTGTAATCCAACCATATCCCTGTA",
            alt="G",
            qual="528",
            filter="PASS",
            info="END=34102668;SVTYPE=DEL;SVLEN=-95;CIGAR=1M95D;CONTIG=ACCATAGTCCTCAAAGCTTTAAATGAATGGGTACCCAGATGCCTGTCCACCTCTTCTCCTATCTCTTGTCCCCTCATCCCTGTGATCCAACCATATCCCTGTGATCCGGCCATATCCCTGTAATCCGGCCATGCTGGACCTTTGCATTTGCCAAGACTTCAGCCCAGAAAATGTTTCCCCAGTAGCCATGGGACCCTTACTTATTTCCGGTGTCTA;CIPOS=0,4;HOMLEN=4;HOMSEQ=ATCC",
            format="GT:FT:GQ:PL:PR:SR",
            genotypes=["0/1:PASS:55:578,0,52:14,0:30,16"],
        )

        alignments = [
            AlignedSegment.from_dict(
                {
                    "name": "MantaDEL:2026:0:0:0:2:0",
                    "flag": "0",
                    "ref_name": "chr1",
                    "ref_pos": "34102641",
                    "map_quality": "60",
                    "cigar": "75S141M",
                    "next_ref_name": "*",
                    "next_ref_pos": "0",
                    "length": "0",
                    "seq": "ACCATAGTCCTCAAAGCTTTAAATGAATGGGTACCCAGATGCCTGTCCACCTCTTCTCCTATCTCTTGTCCCCTCATCCCTGTGATCCAACCATATCCCTGTGATCCGGCCATATCCCTGTAATCCGGCCATGCTGGACCTTTGCATTTGCCAAGACTTCAGCCCAGAAAATGTTTCCCCAGTAGCCATGGGACCCTTACTTATTTCCGGTGTCTA",
                    "qual": "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                    "tags": [
                        "NM:i:2",
                        "ms:i:262",
                        "AS:i:262",
                        "nn:i:0",
                        "tp:A:P",
                        "cm:i:28",
                        "s1:i:173",
                        "s2:i:0",
                        "de:f:0.0142",
                        "rl:i:0",
                    ],
                },
                header=self.header,
            )
        ]

        self.assertFalse(
            check_contig_support(
                variant=variant,
                alignments=alignments,
            )
        )
