from typing import Dict, Iterable

from mantatools.core import Variant


def parse_vcf(stream: Iterable) -> Dict[str, Variant]:
    raise NotImplementedError()
