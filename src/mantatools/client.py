import gzip

import click

from mantatools.parsers import parse_vcf


@click.group()
def client():
    pass


@client.command()
@click.option("--vcf", type=click.Path(exists=True), required=True)
def create_bedpe(vcf: str):
    with gzip.open(vcf, "rt") as stream:
        for variant in parse_vcf(stream).values():
            print(variant.to_bedpe())


def run():
    client()
