import unittest

from mantatools.core import Position, Variant
from mantatools.exceptions import InfoFieldNotFound, GenotypeFieldNotFound


class TestVariant(unittest.TestCase):

    def setUp(self) -> None:
        self.variant = Variant(
            chrom="chr1",
            pos="100",
            id="MyVariant",
            ref="A",
            alt="<DEL>",
            qual="1000",
            filter="PASS",
            info="IMPRECISE;END=200;SVTYPE=DEL",
            format="GT:PR",
            genotypes=["0/1:20,15"],
        )

    def test_get_info(self) -> None:
        self.assertTrue(self.variant.get_info("IMPRECISE"))
        self.assertEqual(self.variant.get_info("END"), "200")
        self.assertEqual(self.variant.get_info("SVTYPE"), "DEL")

    def test_info_field_not_found(self) -> None:
        with self.assertRaises(InfoFieldNotFound):
            self.variant.get_info("NON_EXISTENT_INFO_FIELD")

    def test_get_genotype(self) -> None:
        self.assertEqual(self.variant.get_genotype("GT"), "0/1")
        self.assertEqual(self.variant.get_genotype("PR"), "20,15")

    def test_genotype_field_not_found(self) -> None:
        with self.assertRaises(GenotypeFieldNotFound):
            self.variant.get_genotype("NON_EXISTENT_GENOTYPE_FIELD")


class TestBreakpoints(unittest.TestCase):

    def test_deletion(self) -> None:
        variant = Variant(
            chrom="chr1",
            pos="100",
            id="MantaDEL",
            ref="A",
            alt="<DEL>",
            qual="1000",
            filter="PASS",
            info="END=200;SVTYPE=DEL",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr1", pos=100),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr1", pos=200),
        )

    def test_duplication(self) -> None:
        variant = Variant(
            chrom="chr2",
            pos="200",
            id="MantaDUP",
            ref="C",
            alt="<DUP>",
            qual="1000",
            filter="PASS",
            info="END=300;SVTYPE=DUP",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr2", pos=200),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr2", pos=300),
        )

    def test_inversion(self) -> None:
        variant = Variant(
            chrom="chr3",
            pos="300",
            id="MantaINV",
            ref="G",
            alt="<INV>",
            qual="1000",
            filter="PASS",
            info="END=400;SVTYPE=INV",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr3", pos=300),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr3", pos=400),
        )

    def test_insertion(self) -> None:
        variant = Variant(
            chrom="chr4",
            pos="400",
            id="MantaINS",
            ref="T",
            alt="<INS>",
            qual="1000",
            filter="PASS",
            info="END=400;SVTYPE=DUP",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr4", pos=400),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr4", pos=400),
        )

    def test_breakend(self) -> None:
        variant = Variant(
            chrom="chr5",
            pos="500",
            id="MantaBND:0",
            ref="A",
            alt="]chr6:600]A",
            qual="1000",
            filter="PASS",
            info="SVTYPE=BND;MATEID=MantaBND:1",
            format="GT",
            genotypes=["1/1"],
        )

        variant.mate = Variant(
            chrom="chr6",
            pos="600",
            id="MantaBND:1",
            ref="C",
            alt="C[chr5:500[",
            qual="1000",
            filter="PASS",
            info="SVTYPE=BND;MATEID=MantaBND:0",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr5", pos=500),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr6", pos=600),
        )
