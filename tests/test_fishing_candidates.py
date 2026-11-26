from foam import candidates as cand
from foam import fishing


def paf(rid, qlen, qs, qe, matches, block):
    return "\t".join(map(str, [rid, qlen, qs, qe, "+", "bait", 30000, 0, qe - qs,
                               matches, block, 60])) + "\n"


def test_fish_filter_keeps_only_mostly_mito_reads():
    lines = [
        paf("mito_full", 10000, 0, 10000, 9500, 10000),     # keep
        paf("numt_read", 40000, 5000, 7000, 1900, 2000),    # nuclear read with NUMT: drop
        paf("split", 8000, 0, 4100, 3500, 4100),            # split alignment, union 95%: keep
        paf("split", 8000, 4000, 7600, 3000, 3600),
        paf("short", 500, 0, 500, 500, 500),                # below min_len: drop
        paf("lowid", 5000, 0, 5000, 3000, 5000),            # 60% identity: drop
    ]
    kept = list(fishing.filter_paf(lines, min_frac=0.8, min_aln=1000,
                                   min_identity=0.7, min_len=1000, max_len=0))
    assert kept == ["mito_full", "split"]


def test_max_len():
    kept = list(fishing.filter_paf([paf("r", 50000, 0, 50000, 49000, 50000)], max_len=40000))
    assert kept == []


def test_double_bait():
    assert fishing.double_bait([("b", "", "ACGT")])[0][2] == "ACGTACGT"


def test_flye_info_and_hifiasm_names(tmp_path):
    info = tmp_path / "assembly_info.txt"
    info.write_text("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
                    "contig_1\t15000\t300\tY\tN\t1\t*\t1\n"
                    "contig_2\t30100\t150\tY\tN\t1\t*\t2\n"
                    "contig_3\t15000\t300\tN\tN\t1\t*\t3\n")
    circ = cand.flye_circular(str(info))
    assert circ == {"contig_1": True, "contig_2": True, "contig_3": False}
    assert cand.is_hifiasm_circular("ptg000982c")
    assert cand.is_hifiasm_circular("h1tg000012c")
    assert not cand.is_hifiasm_circular("ptg000982l")


def test_select_lengths_circularity_names():
    recs = [("ptg1c", "", "A" * 15000), ("ptg2c", "", "A" * 30000),
            ("ptg3l", "", "A" * 15000), ("ptg4c", "", "A" * 900000)]
    recs = [("ptg00000%dc" % i if n.endswith("c") else "ptg00000%dl" % i, d, s)
            for i, (n, d, s) in enumerate(recs)]
    keep, table = cand.select(recs, min_len=11000, max_len=75000)
    assert [k[0] for k in keep] == ["ptg000000c", "ptg000001c"]
    reasons = dict((t["contig"], t["reasons"]) for t in table)
    assert reasons["ptg000002l"] == ["not_circular"]
    assert reasons["ptg000003c"] == ["too_long"]
    keep, _ = cand.select(recs, 11000, 75000, names=["ptg000002l"])
    assert [k[0] for k in keep] == ["ptg000002l"]
