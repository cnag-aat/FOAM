import pytest

from conftest import needs_minimap2, rand_seq
from foam import seqio, tandem

P = 6000


def is_rotation(a, b):
    return len(a) == len(b) and a in b + b


@pytest.mark.parametrize("n", [2, 3])
def test_kmer_split_circular(n):
    mono = rand_seq(P, 3)
    contig = seqio.rotate(mono * n, 1234)        # arbitrary linearisation point
    copies, info = tandem.split_copies(contig, True, min_period=2000, use_minimap2=False)
    assert info["period"] == P and info["n_copies"] == n
    for _s, _e, c in copies:
        assert is_rotation(c, mono)


def test_copies_cut_at_homologous_points_with_indel():
    mono = rand_seq(P, 4)
    c2 = mono[:3000] + mono[3001:]               # 1 bp deletion in copy 2
    copies, info = tandem.split_copies(mono + c2, True, min_period=2000, use_minimap2=False)
    assert info["boundaries"] == "anchor"
    assert sorted(e - s for s, e, _ in copies) == [P - 1, P]
    # both copies start at the same biological position
    starts = set(c[:30] for _s, _e, c in copies)
    assert len(starts) == 1


def test_single_copy_untouched():
    s = rand_seq(14000, 5)
    copies, info = tandem.split_copies(s, True, min_period=2000, use_minimap2=False)
    assert info["n_copies"] == 1 and copies[0][2] == s


def test_control_region_repeat_is_not_a_period():
    # short tandem repeat (e.g. in the control region) must not split the genome
    s = rand_seq(14000, 6)
    s = s[:7000] + rand_seq(150, 7) * 8 + s[7000:]
    copies, info = tandem.split_copies(s, True, min_period=5000, use_minimap2=False)
    assert info["n_copies"] == 1


def test_linear_contig_keeps_complete_copies_only():
    mono = rand_seq(P, 8)
    contig = mono[2000:] + mono + mono[:3000]     # partial + full + partial
    copies, info = tandem.split_copies(contig, False, min_period=2000, use_minimap2=False)
    assert copies and all(is_rotation(c, mono) for _s, _e, c in copies)


@needs_minimap2
def test_minimap2_period():
    mono = rand_seq(P, 9)
    period, method = tandem.estimate_period(mono * 2, min_period=2000)
    assert method == "minimap2" and period == P
