import os
import random
import shutil
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)


def rand_seq(n, seed=0):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


@pytest.fixture
def root():
    return ROOT


needs_minimap2 = pytest.mark.skipif(shutil.which("minimap2") is None,
                                    reason="minimap2 not on PATH")
needs_snakemake = pytest.mark.skipif(shutil.which("snakemake") is None,
                                     reason="snakemake not on PATH")
