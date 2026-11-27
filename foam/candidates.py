"""Pick candidate mitochondrial contigs from a flye or hifiasm assembly."""

import os
import re
import shutil
import subprocess

from .fishing import _union_len

HIFIASM_CIRC = re.compile(r"^(?:ptg|atg|utg|h\d+tg)\d+c$")


def flye_circular(info_path):
    """Return {contig: is_circular} from flye's assembly_info.txt."""
    circ = {}
    with open(info_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            circ[c[0]] = c[3].strip().upper() == "Y"
    return circ


def is_hifiasm_circular(name):
    """hifiasm marks circular contigs with a trailing 'c' (ptg000012c)."""
    return bool(HIFIASM_CIRC.match(name))


def bait_coverage(contigs_fa, bait_fa, minimap2="minimap2", threads=4, preset="asm20"):
    """Fraction of each contig covered by alignments to the bait."""
    exe = shutil.which(minimap2) if os.path.sep not in minimap2 else minimap2
    if not exe:
        raise RuntimeError("minimap2 not found; needed for bait screening")
    p = subprocess.Popen([exe, "-x", preset, "-t", str(threads), bait_fa, contigs_fa],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        raise RuntimeError("minimap2 failed: %s" % err.decode(errors="replace"))
    ivs, lens = {}, {}
    for line in out.decode().splitlines():
        c = line.split("\t")
        if len(c) < 12:
            continue
        lens[c[0]] = int(c[1])
        ivs.setdefault(c[0], []).append((int(c[2]), int(c[3])))
    return dict((k, _union_len(v) / float(lens[k])) for k, v in ivs.items())


def select(records, min_len, max_len, circular_map=None, require_circular=True,
           names=None, bait_cov=None, min_bait_cov=0.3):
    """Return (selected records, per-contig decision table)."""
    keep, table = [], []
    wanted = set(names) if names else None
    for name, _desc, seq in records:
        L = len(seq)
        circ = circular_map.get(name, False) if circular_map is not None else is_hifiasm_circular(name)
        reasons = []
        if wanted is not None:
            if name not in wanted:
                reasons.append("not_requested")
        else:
            if L < min_len:
                reasons.append("too_short")
            if max_len and L > max_len:
                reasons.append("too_long")
            if require_circular and not circ:
                reasons.append("not_circular")
            if bait_cov is not None and bait_cov.get(name, 0.0) < min_bait_cov:
                reasons.append("no_bait_hit")
        table.append({"contig": name, "length": L, "circular": circ,
                      "bait_cov": None if bait_cov is None else round(bait_cov.get(name, 0.0), 3),
                      "selected": not reasons, "reasons": reasons})
        if not reasons:
            keep.append((name, "circular" if circ else "linear", seq))
    return keep, table
