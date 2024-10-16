import gzip

import click

from collections import defaultdict
from typing import Dict, List, Optional

from pysam import AlignedSegment, AlignmentFile

from svtoolbox.exceptions import InfoFieldNotFound, SVToolBoxException
from svtoolbox.parser import parse_vcf


@click.group()
def client():
    pass


@client.command()
@click.option("--vcf", type=click.Path(exists=True), required=True)
@click.option("--include_fields", type=str, required=False)
def create_bedpe(vcf: str, include_fields: Optional[str] = None) -> None:
    with gzip.open(vcf, "rt") as stream:
        for variant in parse_vcf(stream).values():
            print(
                variant.to_bedpe(
                    include_fields=include_fields.split(",") if include_fields else None
                )
            )


@client.command()
@click.option("--vcf", type=click.Path(exists=True), required=True)
def create_contigs_fastq(vcf: str) -> None:
    with gzip.open(vcf, "rt") as stream:
        for variant in parse_vcf(stream).values():
            try:
                contig = str(variant.get_info("CONTIG"))
                print("\n".join([f"@{variant.id}", contig, "+", "I" * len(contig)]))
            except InfoFieldNotFound:
                pass


def run():
    client()
