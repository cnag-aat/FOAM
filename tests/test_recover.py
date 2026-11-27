from conftest import rand_seq
from foam import annotation as ann
from foam import recover, seqio
from foam.gencode import GeneticCode


def stopless(n, seed, code=5):
    """Random coding sequence of n codons without stops (start ATG)."""
    gc = GeneticCode(code)
    s = rand_seq(n * 3 * 3, seed)
    cods = [s[i:i + 3] for i in range(0, len(s), 3)]
    cods = [c for c in cods if c not in gc.stops and c not in gc.starts][:n - 1]
    return "ATG" + "".join(cods)


def test_genetic_codes():
    gc5 = GeneticCode(5)
    assert gc5.codon["AGA"] == "S" and gc5.codon["TGA"] == "W" and gc5.codon["ATA"] == "M"
    assert {"ATT", "ATA", "ATG", "TTG", "GTG"} <= gc5.starts
    assert GeneticCode(2).codon["AGA"] == "*"


def _genome_with_orf(strand):
    orf = stopless(100, 11) + "TAA"
    left, right = "TAG" * 20 + rand_seq(37, 1), rand_seq(40, 2) + "TAA" * 20
    gene = orf if strand > 0 else seqio.revcomp(orf)
    if strand < 0:
        left, right = seqio.revcomp(right), seqio.revcomp(left)
    seq = left + gene + right
    return seq, len(left), len(left) + len(orf)


def test_orf_from_hit_both_strands():
    for strand in (1, -1):
        seq, s, e = _genome_with_orf(strand)
        # hit covering codons 10..80 of the protein
        if strand > 0:
            hit = {"start": s + 30, "end": s + 240, "strand": 1, "sstart": 11}
        else:
            hit = {"start": e - 1 - 30, "end": e - 1 - 240, "strand": -1, "sstart": 11}
        f, prot = recover.orf_from_hit(seq, hit, 5)
        assert (f.start, f.end, f.strand) == (s, e, strand)
        assert prot.startswith("M") and prot.endswith("*") and prot.count("*") == 1


def test_trim_rrna_to_neighbours():
    L = 5000
    feats = [ann.Feature("trnL1", "tRNA", 1000, 1060, -1),
             ann.Feature("trnV", "tRNA", 2040, 2100, -1)]
    r = ann.Feature("rrnL", "rRNA", 1040, 2080, -1, source="ignored.mitos")
    assert recover.trim_to_neighbours(r, feats, L)
    assert (r.start, r.end) == (1060, 2040)



def test_trim_rrna_order_independent_and_contained():
    # xgXerMont: cox1 | trnV | rrnL; the untrimmed rrnL starts inside cox1 and
    # runs over all of trnV. Trimming to cox1 alone left trnV inside rrnL.
    L = 14024
    for order in (1, -1):
        feats = [ann.Feature("cox1", "gene", 0, 1530, 1),
                 ann.Feature("trnV", "tRNA", 1526, 1587, 1),
                 ann.Feature("trnL1", "tRNA", 2598, 2659, 1)][::order]
        r = ann.Feature("rrnL", "rRNA", 1506, 2598, 1, source="ignored.mitos")
        assert recover.trim_to_neighbours(r, feats, L)
        assert (r.start, r.end) == (1587, 2598)
    # a tRNA wholly inside, near the 3' end: cut on that side
    feats = [ann.Feature("trnX", "tRNA", 1900, 1960, 1)]
    r = ann.Feature("rrnL", "rRNA", 1000, 2000, 1)
    assert recover.trim_to_neighbours(r, feats, L) and (r.start, r.end) == (1000, 1900)
    # a gene in the middle: too much would be lost
    feats = [ann.Feature("nad1", "gene", 1400, 1600, 1)]
    r = ann.Feature("rrnL", "rRNA", 1000, 2000, 1)
    assert not recover.trim_to_neighbours(r, feats, L)


def test_recover_rrna_from_ignored():
    L = 5000
    feats = [ann.Feature("trnL1", "tRNA", 1000, 1060, -1)]
    cands = [ann.Feature("rrnL", "rRNA", 1040, 2080, -1, score=0.0, method="mitfi",
                         source="ignored.mitos")]
    log = recover.recover_rrnas(feats, cands, ["rrnL", "trnL1"], L)
    assert log and [f.name for f in feats] == ["trnL1", "rrnL"]


def test_fix_5prime_overlap_moves_start():
    # + strand protein starting at 100; tRNA covers 70..130; in-frame ATG at 130
    body = stopless(60, 12)
    seq = rand_seq(100, 13) + "ATA" + "CTT" * 9 + "ATG" + body[3:] + "TAA" + rand_seq(200, 14)
    p = ann.Feature("nad2", "gene", 100, 100 + 3 + 27 + len(body) + 3, 1)
    t = ann.Feature("trnI", "tRNA", 70, 130, 1)
    log = recover.fix_5prime_overlaps([p, t], seq, 5, len(seq))
    assert log and p.start == 130
