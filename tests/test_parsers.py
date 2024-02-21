import unittest

from mantatools.core import Position
from mantatools.parsers import parse_vcf


class TestVcfParser(unittest.TestCase):

    def setUp(self) -> None:

        vcf_lines = [
            "##fileformat=VCFv4.1",
            '##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">',
            '##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999">',
            "\t".join(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "NORMAL",
                    "TUMOR",
                ]
            ),
            "\t".join(
                [
                    "chr1",
                    "100",
                    "MantaDEL",
                    "A",
                    "<DEL>",
                    ".",
                    "PASS",
                    "END=200;SVTYPE=DEL",
                    "PR:SR",
                    "15,0:30,0",
                    "30,5:60,20",
                ]
            ),
            "\t".join(
                [
                    "chr2",
                    "200",
                    "MantaBND:0",
                    "A",
                    "[chr4:400[A",
                    ".",
                    "PASS",
                    "SVTYPE=BND;MATEID=MantaBND:1",
                    "PR:SR",
                    "25,0:15,0",
                    "20,10:50,30",
                ]
            ),
            "\t".join(
                [
                    "chr3",
                    "300",
                    "MantaDUP",
                    "G",
                    "<DUP>",
                    ".",
                    "PASS",
                    "END=400;SVTYPE=BND",
                    "PR:SR",
                    "45,0:30,0",
                    "25,15:45,45",
                ]
            ),
            "\t".join(
                [
                    "chr4",
                    "400",
                    "MantaBND:1",
                    "G",
                    "[chr2:200[G",
                    ".",
                    "PASS",
                    "SVTYPE=BND;MATEID=MantaBND:0",
                    "PR:SR",
                    "25,0:15,0",
                    "20,10:50,30",
                ]
            ),
        ]

        self.variants = parse_vcf(vcf_lines)

    def test_variant_id(self) -> None:
        self.assertEqual(self.variants[0].id, "MantaDEL")
        self.assertEqual(self.variants[1].id, "MantaBND:0")
        self.assertEqual(self.variants[2].id, "MantaDUP")
        self.assertEqual(self.variants[3].id, "MantaBND:1")

    def test_variant_start(self) -> None:
        self.assertEqual(self.variants[0].start, Position(chrom="chr1", pos=100))
        self.assertEqual(self.variants[1].start, Position(chrom="chr2", pos=200))
        self.assertEqual(self.variants[2].start, Position(chrom="chr3", pos=300))
        self.assertEqual(self.variants[3].start, Position(chrom="chr4", pos=400))

    def test_variant_end(self) -> None:
        self.assertEqual(self.variants[0].end, Position(chrom="chr1", pos=200))
        self.assertEqual(self.variants[1].end, Position(chrom="chr4", pos=400))
        self.assertEqual(self.variants[2].end, Position(chrom="chr3", pos=400))
        self.assertEqual(self.variants[3].end, Position(chrom="chr2", pos=200))
