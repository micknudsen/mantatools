import unittest

from svtoolbox.core import BedPE, Interval, Position, Variant
from svtoolbox.exceptions import (
    FieldNotFound,
    InfoFieldNotFound,
    GenotypeFieldNotFound,
    MissingMate,
)


class TestInterval(unittest.TestCase):

    def test_overlaps(self) -> None:
        interval_1 = Interval(chrom="chr1", left=100, right=200)
        interval_2 = Interval(chrom="chr1", left=150, right=250)
        interval_3 = Interval(chrom="chr1", left=250, right=300)
        interval_4 = Interval(chrom="chr2", left=100, right=200)

        self.assertTrue(interval_1.overlaps(interval_2))
        self.assertFalse(interval_1.overlaps(interval_3))
        self.assertFalse(interval_1.overlaps(interval_4))
        self.assertTrue(interval_2.overlaps(interval_3))
        self.assertFalse(interval_2.overlaps(interval_4))
        self.assertFalse(interval_3.overlaps(interval_4))


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
            info="IMPRECISE;END=200;SVTYPE=DEL;CIPOS=-10,5;CIEND=-15,20",
            format="GT:PR",
            genotypes={"SAMPLE": "0/1:20,15"},
        )

    def test_get_info(self) -> None:
        self.assertTrue(self.variant.get_info("IMPRECISE"))
        self.assertEqual(self.variant.get_info("END"), "200")
        self.assertEqual(self.variant.get_info("SVTYPE"), "DEL")
        self.assertEqual(self.variant.get_info("CIPOS"), "-10,5")
        self.assertEqual(self.variant.get_info("CIEND"), "-15,20")

    def test_set_info(self) -> None:
        self.variant.set_info(key="SOME_NEW_FIELD", value="SOME_NEW_VALUE")
        self.assertEqual(self.variant.get_info("SOME_NEW_FIELD"), "SOME_NEW_VALUE")

    def test_info_field_not_found(self) -> None:
        with self.assertRaises(InfoFieldNotFound):
            self.variant.get_info("NON_EXISTENT_INFO_FIELD")

    def test_get_genotype(self) -> None:
        self.assertEqual(self.variant.get_genotype(sample="SAMPLE", key="GT"), "0/1")
        self.assertEqual(self.variant.get_genotype(sample="SAMPLE", key="PR"), "20,15")

    def test_genotype_field_not_found(self) -> None:
        with self.assertRaises(GenotypeFieldNotFound):
            self.variant.get_genotype(
                sample="SAMPLE", key="NON_EXISTENT_GENOTYPE_FIELD"
            )

    def test_to_bedpe(self) -> None:
        self.assertEqual(
            self.variant.to_bedpe(),
            BedPE(
                chrom_1="chr1",
                start_1=89,
                end_1=105,
                chrom_2="chr1",
                start_2=184,
                end_2=220,
                name="MyVariant",
                score="1000",
                strand_1=None,
                strand_2=None,
            ),
        )

    def test_to_bedpe_with_extra_fields(self) -> None:
        self.assertEqual(
            self.variant.to_bedpe(
                include_fields=[
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                ]
            ),
            BedPE(
                chrom_1="chr1",
                start_1=89,
                end_1=105,
                chrom_2="chr1",
                start_2=184,
                end_2=220,
                name="MyVariant",
                score="1000",
                strand_1=None,
                strand_2=None,
                fields={
                    "REF": "A",
                    "ALT": "<DEL>",
                    "QUAL": "1000",
                    "FILTER": "PASS",
                },
            ),
        )

    def test_to_bedpe_field_not_found(self) -> None:
        with self.assertRaises(FieldNotFound):
            self.variant.to_bedpe(include_fields=["NON_EXISTENT_FIELD"])

    def test_str_method(self) -> None:
        self.assertEqual(
            self.variant.__str__(),
            "chr1\t100\tMyVariant\tA\t<DEL>\t1000\tPASS\tIMPRECISE;END=200;SVTYPE=DEL;CIPOS=-10,5;CIEND=-15,20\tGT:PR\t0/1:20,15",
        )

    def test_str_method_after_insertion_new_key_value_info(self) -> None:
        self.variant.set_info(key="NEW_FIELD", value="NEW_VALUE")
        self.assertEqual(
            self.variant.__str__(),
            "chr1\t100\tMyVariant\tA\t<DEL>\t1000\tPASS\tIMPRECISE;END=200;SVTYPE=DEL;CIPOS=-10,5;CIEND=-15,20;NEW_FIELD=NEW_VALUE\tGT:PR\t0/1:20,15",
        )

    def test_str_method_after_insertion_new_flag_info(self) -> None:
        self.variant.set_info(key="NEW_FLAG", value=True)
        self.assertEqual(
            self.variant.__str__(),
            "chr1\t100\tMyVariant\tA\t<DEL>\t1000\tPASS\tIMPRECISE;END=200;SVTYPE=DEL;CIPOS=-10,5;CIEND=-15,20;NEW_FLAG\tGT:PR\t0/1:20,15",
        )


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
            genotypes={"SAMPLE": "1/1"},
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
            genotypes={"SAMPLE": "1/1"},
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
            genotypes={"SAMPLE": "1/1"},
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
            genotypes={"SAMPLE": "1/1"},
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr4", pos=400),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr4", pos=400),
        )

    def test_breakend_manta_style(self) -> None:
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
            genotypes={"SAMPLE": "1/1"},
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
            genotypes={"SAMPLE": "1/1"},
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr6", pos=600),
        )

    def test_breakend_delly_style(self) -> None:
        variant = Variant(
            chrom="chr7",
            pos="2000",
            id="BND000012345",
            ref="A",
            alt="A]chr8:3000]",
            qual="1000",
            filter="PASS",
            info="SVTYPE=BND;CHR2=chr8;POS2=3000",
            format="GT",
            genotypes={"SAMPLE": "1/1"},
        )

        self.assertEqual(
            variant.start,
            Position(chrom="chr7", pos=2000),
        )

        self.assertEqual(
            variant.end,
            Position(chrom="chr8", pos=3000),
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
            genotypes={"SAMPLE": "1/1"},
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr7", left=690, right=705),
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
            genotypes={"SAMPLE": "1/1"},
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
            genotypes={"SAMPLE": "1/1"},
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
            genotypes={"SAMPLE": "1/1"},
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr10", left=1000, right=1000),
        )

        self.assertEqual(
            variant.ci_end,
            Interval(chrom="chr10", left=1100, right=1100),
        )

    def test_ci_start_and_ci_end_breakend_manta_style(self) -> None:
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
            genotypes={"SAMPLE": "1/1"},
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
            genotypes={"SAMPLE": "1/1"},
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr5", left=500, right=500),
        )
        self.assertEqual(
            variant.ci_end,
            Interval(chrom="chr6", left=600, right=600),
        )

    def test_ci_start_and_ci_end_breakend_delly_style(self) -> None:
        variant = Variant(
            chrom="chr7",
            pos="700",
            id="BND000012345",
            ref="A",
            alt="]chr8:800]A",
            qual="1000",
            filter="PASS",
            info="SVTYPE=BND;CHR2=chr8;POS2=800;CIPOS=-10,10;CIEND=-20,20",
            format="GT",
            genotypes={"SAMPLE": "1/1"},
        )

        self.assertEqual(
            variant.ci_start,
            Interval(chrom="chr7", left=690, right=710),
        )
        self.assertEqual(
            variant.ci_end,
            Interval(chrom="chr8", left=780, right=820),
        )


class TestBedPE(unittest.TestCase):

    def test_str_method_with_optional_fields(self) -> None:
        bedpe = BedPE(
            chrom_1="chr1",
            start_1=100,
            end_1=200,
            chrom_2="chr3",
            start_2=400,
            end_2=500,
            name="MyVariant",
            score="1000",
            strand_1="+",
            strand_2="-",
        )

        self.assertEqual(
            bedpe.__str__(),
            "chr1\t100\t200\tchr3\t400\t500\tMyVariant\t1000\t+\t-",
        )

    def test_str_method_without_optional_fields(self) -> None:
        bedpe = BedPE(
            chrom_1="chr1",
            start_1=100,
            end_1=200,
            chrom_2="chr3",
            start_2=400,
            end_2=500,
            name="MyVariant",
            score=None,
            strand_1=None,
            strand_2=None,
        )

        self.assertEqual(
            bedpe.__str__(),
            "chr1\t100\t200\tchr3\t400\t500\tMyVariant\t.\t.\t.",
        )

    def test_str_method_with_extra_fields(self) -> None:
        bedpe = BedPE(
            chrom_1="chr1",
            start_1=100,
            end_1=200,
            chrom_2="chr3",
            start_2=400,
            end_2=500,
            name="MyVariant",
            score="1000",
            strand_1="+",
            strand_2="-",
            fields={
                "REF": "A",
                "ALT": "<DEL>",
                "QUAL": "1000",
                "FILTER": "PASS",
            },
        )

        self.assertEqual(
            bedpe.__str__(),
            "chr1\t100\t200\tchr3\t400\t500\tMyVariant\t1000\t+\t-\tREF=A;ALT=<DEL>;QUAL=1000;FILTER=PASS",
        )
