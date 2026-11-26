import random

from conftest import rand_seq
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


def test_kmer_select_keeps_mito_rejects_numt_reads(tmp_path):
    path = simulate(tmp_path)
    keep, rep = kmers.select([path], k=21, mod=4)
    kept = set(keep)
    mito = [r for r in kept if r.startswith("mito_")]
    numt = [r for r in kept if r.startswith("numt_")]
    assert len(mito) >= 0.97 * 300
    assert numt == []
    # threshold sits between nuclear (NUMT locus) coverage and the mito peak
    assert 4 < rep["abundance_threshold"] < rep["coverage_peak"]


def test_kmer_select_independent_of_read_accuracy(tmp_path):
    for err in (0.005, 0.06):
        keep, rep = kmers.select([simulate(tmp_path, error=err)], k=21, mod=4)
        assert sum(r.startswith("mito_") for r in keep) >= 0.95 * 300
        assert not any(r.startswith("numt_") for r in keep)


def test_fastq_reading(tmp_path):
    p = tmp_path / "r.fq"
    p.write_text("@a x\nACGT\n+\nIIII\n@b\nGGCC\n+\nIIII\n")
    assert list(kmers.iter_reads(str(p))) == [("a", "ACGT"), ("b", "GGCC")]


def test_no_mito_reports_error(tmp_path):
    p = tmp_path / "r.fa"
    seqio.write_fasta([("r%d" % i, rand_seq(2000, i)) for i in range(20)], str(p))
    keep, rep = kmers.select([str(p)])
    assert keep == [] and "error" in rep
