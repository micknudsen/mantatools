import unittest

from mantatools.core import Variant
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
            chrom="chr9",
            pos="22496528",
            id="MantaDEL",
            ref="T",
            alt="<DEL>",
            qual="999",
            filter="PASS",
            info="END=22504345;SVTYPE=DEL",
            format="GT",
            genotypes=["1/1"],
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr9", pos=22496528),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr9", pos=22504345),
        )
