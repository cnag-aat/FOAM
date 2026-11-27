import collections
import os
import random

from conftest import needs_fastk, rand_seq
from foam import kmers, seqio


def mutate(seq, rate, rnd):
    out = []
    for c in seq:
        r = rnd.random()
        if r < rate / 3:
            continue                                # deletion
        if r < 2 * rate / 3:
            out.append(rnd.choice("ACGT"))          # substitution
        elif r < rate:
            out.extend([c, rnd.choice("ACGT")])     # insertion
        else:
            out.append(c)
    return "".join(out)


def simulate(tmp_path, error=0.02, mito_reads=300, numt_loci=15, nuc_cov=4):
    rnd = random.Random(7)
    mito = rand_seq(6000, 1)
    reads = []
    for i in range(mito_reads):
        start, ln = rnd.randrange(6000), rnd.randint(1500, 4000)
        s = (mito + mito)[start:start + ln]
        if rnd.random() < 0.5:
            s = seqio.revcomp(s)
        reads.append(("mito_%d" % i, mutate(s, error, rnd)))
    for j in range(numt_loci):
        # nuclear locus with a 1.2 kb NUMT, 97% identical to the mitogenome
        a = rnd.randrange(6000)
        numt = mutate((mito + mito)[a:a + 1200], 0.03, rnd)
        locus = rand_seq(3000, 100 + j) + numt + rand_seq(3000, 200 + j)
        for c in range(nuc_cov):
            st = rnd.randrange(0, 1500)
            reads.append(("numt_%d_%d" % (j, c), mutate(locus[st:st + 5500], error, rnd)))
    path = tmp_path / "cand.fa"
    seqio.write_fasta(reads, str(path))
    return str(path)


import pytest

BACKENDS = ["python", pytest.param("fastk", marks=needs_fastk)]


@pytest.mark.parametrize("backend", BACKENDS)
def test_kmer_select_keeps_mito_rejects_numt_reads(tmp_path, backend):
    path = simulate(tmp_path)
    keep, rep = kmers.select([path], k=21, mod=4, backend=backend, tmpdir=str(tmp_path))
    assert rep["backend"] == backend
    kept = set(keep)
    mito = [r for r in kept if r.startswith("mito_")]
    numt = [r for r in kept if r.startswith("numt_")]
    assert len(mito) >= 0.97 * 300
    assert numt == []
    # threshold sits between nuclear (NUMT locus) coverage and the mito peak
    assert 4 < rep["abundance_threshold"] < rep["coverage_peak"]


@pytest.mark.parametrize("backend", BACKENDS)
def test_kmer_select_independent_of_read_accuracy(tmp_path, backend):
    for err in (0.005, 0.06):
        keep, rep = kmers.select([simulate(tmp_path, error=err)], k=21, mod=4,
                                 backend=backend, tmpdir=str(tmp_path))
        assert sum(r.startswith("mito_") for r in keep) >= 0.95 * 300
        assert not any(r.startswith("numt_") for r in keep)


def test_fastq_reading(tmp_path):
    p = tmp_path / "r.fq"
    p.write_text("@a x\nACGT\n+\nIIII\n@b\nGGCC\n+\nIIII\n")
    assert list(kmers.iter_reads(str(p))) == [("a", "ACGT"), ("b", "GGCC")]


@pytest.mark.parametrize("backend", BACKENDS)
def test_no_mito_reports_error(tmp_path, backend):
    p = tmp_path / "r.fa"
    seqio.write_fasta([("r%d" % i, rand_seq(2000, i)) for i in range(20)], str(p))
    keep, rep = kmers.select([str(p)], backend=backend, tmpdir=str(tmp_path))
    assert keep == [] and "error" in rep


def test_peak_from_hist_matches_counts():
    rnd = random.Random(3)
    counts = {i: rnd.choice([1, 1, 2, 3, 150, 160, 170, 180]) for i in range(5000)}
    h = collections.Counter(counts.values())
    assert kmers.peak_from_hist(h) == kmers.coverage_peak(counts)
    assert kmers.histogram_from_hist(h) == kmers.histogram(counts)


@needs_fastk
def test_fastk_read_order_with_short_and_empty_reads(tmp_path):
    """Profiles are matched to reads by ordinal; FastK drops empty FASTA
    records, so FOAM must keep names and profiles aligned itself."""
    a, b = tmp_path / "a.fa", tmp_path / "b.fq"
    a.write_text(">e1\n>s1\nACGT\n>long1\n%s\n>n1\n%s\n" % ("ACGTTGCA" * 20, "N" * 40))
    b.write_text("@long2\n%s\n+\n%s\n@e2\n\n+\n\n" % ("ACGTTGCA" * 5, "I" * 40))
    fk = kmers.FastK([str(a), str(b)], k=21, tmpdir=str(tmp_path), threads=4)
    try:
        assert fk.short_ids == ["e1", "s1", "e2"]
        rows = list(fk.scores(1, workers=2))
    finally:
        fk.close()
    assert [(r, n) for r, n, _ in rows] == [("long1", 140), ("n1", 0), ("long2", 20)]
    assert not [f for f in os.listdir(str(tmp_path)) if f.startswith("foam_fastk")]
