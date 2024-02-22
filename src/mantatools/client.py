import gzip

import click

from mantatools.exceptions import InfoFieldNotFound
from mantatools.parser import parse_vcf


@click.group()
def client():
    pass


@client.command()
@click.option("--vcf", type=click.Path(exists=True), required=True)
def create_bedpe(vcf: str) -> None:
    with gzip.open(vcf, "rt") as stream:
        for variant in parse_vcf(stream).values():
            print(variant.to_bedpe())


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
