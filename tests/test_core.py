import unittest

from mantatools.core import Position, Variant
from mantatools.exceptions import InfoFieldNotFound, GenotypeFieldNotFound, MissingMate


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

        self.assertEqual(
            variant.start,
            Position(chrom="chr5", pos=500),
        )

        with self.assertRaises(MissingMate):
            variant.end

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
            variant.end,
            Position(chrom="chr6", pos=600),
        )


class TestConfidenceIntervals(unittest.TestCase):

    def test_ci_start_and_ci_end(self) -> None:
        variant = Variant(
            chrom="chr7",
            pos="700",
            id="MantaDEL",
            ref="A",
            alt="<DEL>",
            qual="1000",
            filter="PASS",
            info="END=800;SVTYPE=DEL;CIPOS=-10,5;CIEND=-15,20",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr7", left=790, right=805),
        )

        self.assertEqual(
            variant.ci_end,
            Interval(chrom="chr7", left=785, right=820),
        )

    def test_ci_start_only(self) -> None:
        variant = Variant(
            chrom="chr8",
            pos="800",
            id="MantaDEL",
            ref="A",
            alt="<DEL>",
            qual="1000",
            filter="PASS",
            info="END=900;SVTYPE=DEL;CIPOS=-10,5",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr8", left=790, right=805),
        )

        self.assertEqual(
            variant.ci_end,
            Interval(chrom="chr8", left=900, right=900),
        )

    def test_ci_end_only(self) -> None:
        variant = Variant(
            chrom="chr9",
            pos="900",
            id="MantaDEL",
            ref="A",
            alt="<DEL>",
            qual="1000",
            filter="PASS",
            info="END=1000;SVTYPE=DEL;CIEND=-15,20",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr9", left=900, right=900),
        )

        self.assertEqual(
            variant.ci_end,
            Interval(chrom="chr9", left=985, right=1020),
        )

    def test_no_ci_pos_or_ci_end(self) -> None:
        variant = Variant(
            chrom="chr10",
            pos="1000",
            id="MantaDEL",
            ref="A",
            alt="<DEL>",
            qual="1000",
            filter="PASS",
            info="END=1100;SVTYPE=DEL",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr10", left=1000, right=1000),
        )

        self.assertEqual(
            variant.ci_end,
            Interval(chrom="chr10", left=1100, right=1100),
        )
