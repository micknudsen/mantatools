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

    def test_check_contig_support_imprecise_variant(self) -> None:
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
