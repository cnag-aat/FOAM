import random

import pytest

from foam import annotation as ann
from foam import mitos as mt
from foam import seqio


def rand_seq(n, seed=0):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def feat_seq(seq, f):
    s = (seq + seq)[f.start:f.end]
    return s if f.strand > 0 else seqio.revcomp(s)


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("start", [0, 100, 950, 990])  # 990 wraps the origin
def test_orient_puts_anchor_5prime_at_position_1(strand, start):
    L = 1000
    seq = rand_seq(L, 1)
    cox1 = ann.normalize(ann.Feature("cox1", "gene", start, start + 40, strand), L)
    other = ann.normalize(ann.Feature("nad1", "gene", 500, 560, -strand), L)
    before = {f.name: feat_seq(seq, f) for f in (cox1, other)}
    oseq, ofeats, info = mt.orient(seq, [cox1, other], ["trnF", "cox1"])
    assert info["anchor"] == "cox1"
    o = dict((f.name, f) for f in ofeats)
    assert o["cox1"].start == 0 and o["cox1"].strand == 1
    # every feature still covers exactly the same biological sequence
    for name, s in before.items():
        assert feat_seq(oseq, o[name]) == s


def test_vertebrate_preset_uses_trnF_others_cox1():
    assert ann.taxon_preset("vertebrata")["anchor"][0] == "trnF"
    assert ann.taxon_preset("Vertebrata")["genetic_code"] == 2
    for t in ann.TAXON_PRESETS:
        if t != "vertebrata":
            assert ann.TAXON_PRESETS[t]["anchor"] == ["cox1"]
    assert ann.taxon_preset("mollusca")["genetic_code"] == 5
    with pytest.raises(ValueError):
        ann.taxon_preset("fungi")


def test_anchor_uses_5prime_part_of_split_gene():
    a = ann.Feature("cox1", "gene", 100, 300, 1, part=0)
    b = ann.Feature("cox1", "gene", 305, 600, 1, part=1)
    assert ann.choose_anchor([b, a], ["cox1"]) is a


def test_mitos_wrapping_feature_parsed(tmp_path):
    p = tmp_path / "result.mitos"
    # stop < start -> wraps the origin
    p.write_text("acc\tgene\tcox1\tmitos\t990\t29\t1\t5.0\t.\t-\t-\t.\t.\t.\n")
    f = ann.parse_mitos_file(str(p), 1000)[0]
    assert (f.start, f.end, f.length) == (990, 1030, 40)


def test_leu_ser_relabel():
    fs = [ann.Feature("trnL1", "tRNA", 0, 70, 1, anticodon="TAA"),
          ann.Feature("trnS2", "tRNA", 100, 170, 1, anticodon="GCT"),
          ann.Feature("trnL2", "tRNA", 200, 270, 1, anticodon="TAG")]
    ann.relabel_leu_ser(fs)
    assert [f.name for f in fs] == ["trnL2", "trnS1", "trnL1"]


def test_recover_missing_respects_overlap_and_evalue():
    L = 2000
    feats = [ann.Feature("cox1", "gene", 0, 500, 1)]
    cands = [ann.Feature("trnK", "tRNA", 450, 520, 1, score=1e-5, method="mitfi", source="raw:x"),
             ann.Feature("trnK", "tRNA", 800, 870, 1, score=0.5, method="mitfi", source="raw:x"),
             ann.Feature("trnK", "tRNA", 900, 970, 1, score=1e-4, method="mitfi", source="raw:y")]
    log = ann.recover_missing(feats, cands, ["cox1", "trnK"], L)
    k = [f for f in feats if f.name == "trnK"]
    assert len(k) == 1 and k[0].start == 900 and k[0].source == "recovered:raw:y"
    assert log and "trnK" in log[0]


def test_parse_mitfi_offset(tmp_path):
    d = tmp_path / "m" / "mitfi-global"
    d.mkdir(parents=True)
    head = ">acc|5|11|80|1|1|69|50.0|1e-10|0|40|30|TTC|F|mt-tRNA-Phe-full|1|1|g"
    (d / "sequence.fas-500_tRNAout.nc").write_text(head + "\nACGT\n....\n")
    h = ann.parse_mitfi_dir(str(tmp_path / "m"), 1000)[0]
    assert h.name == "trnF" and h.start == 510 and h.end == 580 and h.score == 1e-10


def test_qc_and_gff(tmp_path):
    L = 3000
    feats = [ann.Feature("cox1", "gene", 0, 1500, 1),
             ann.Feature("trnF", "tRNA", 1600, 1670, -1, anticodon="GAA"),
             ann.Feature("rrnL", "rRNA", 1700, 2900, 1)]
    q = ann.qc(feats, ann.expected_genes("metazoa"), L)
    assert q["n_present"] == 3 and "cox2" in q["missing"]
    out = tmp_path / "x.gff3"
    ann.write_gff3(feats, "sp_MT", L, str(out), "spMT", 5)
    txt = out.read_text()
    assert "Is_circular=true" in txt and "transl_table=5" in txt
    assert "product=tRNA-Phe" in txt and "\tCDS\t1\t1500\t" in txt


def test_bed_splits_origin_crossing_feature(tmp_path):
    L = 14024
    feats = [ann.Feature("cox1", "gene", 0, 1530, 1),
             ann.Feature("trnK", "tRNA", 13970, 14031, 1, anticodon="ttt")]
    p = str(tmp_path / "x.bed")
    ann.write_bed(feats, "mt", p, L)
    rows = [l.split("\t") for l in open(p).read().splitlines()]
    assert [(r[1], r[2], r[3]) for r in rows] == [
        ("0", "7", "trnK(ttt)"), ("0", "1530", "cox1"), ("13970", "14024", "trnK(ttt)")]
    assert all(int(r[2]) <= L for r in rows)
