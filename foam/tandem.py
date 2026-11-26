"""Detect and split tandem copies of a mitogenome inside one contig.

Assemblers (hifiasm in particular, flye with --meta sometimes) can output a
circular contig that is really 2-3 head-to-tail copies of the mitogenome.

1. *Period*: self-align the contig with ``minimap2 -x asm20 -DP``; same-strand
   off-diagonal hits have a shift (target start - query start) that is a
   multiple of the repeat unit.  The smallest shift >= ``min_period`` that
   is supported by a long alignment is the period.  Without minimap2 we fall
   back to the most frequent distance between repeated 25-mers.
2. *Copy boundaries*: an anchor k-mer that occurs once per copy (spaced by
   ~period) gives homologous cut points, so every copy starts at the same
   biological position even if copies differ by small indels.
3. Copies are returned as sequences; for a circular contig the last copy
   wraps around the end.
"""

import collections
import os
import shutil
import subprocess
import tempfile

from . import seqio

K = 25


def _which(prog):
    return shutil.which(prog)


def minimap2_shifts(seq, minimap2="minimap2", threads=1, min_block=1000):
    """Return [(shift, aligned_len)] from a minimap2 self-alignment."""
    exe = _which(minimap2) if not os.path.sep in minimap2 else minimap2
    if not exe:
        return None
    tmpd = tempfile.mkdtemp(prefix="foam_self_")
    try:
        fa = os.path.join(tmpd, "c.fa")
        seqio.write_fasta([("c", seq)], fa)
        p = subprocess.Popen([exe, "-x", "asm20", "-DP", "-t", str(threads), fa, fa],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode != 0:
            raise RuntimeError("minimap2 failed: %s" % err.decode(errors="replace"))
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    shifts = []
    for line in out.decode().splitlines():
        c = line.split("\t")
        if len(c) < 12 or c[4] != "+":
            continue
        qs, ts, block = int(c[2]), int(c[7]), int(c[10])
        if block >= min_block and ts > qs:
            shifts.append((ts - qs, block))
    return shifts


def kmer_shifts(seq, k=K, max_occ=6):
    """Distances between consecutive occurrences of repeated k-mers."""
    pos = collections.defaultdict(list)
    for i in range(len(seq) - k + 1):
        pos[seq[i:i + k]].append(i)
    counts = collections.Counter()
    for plist in pos.values():
        if 1 < len(plist) <= max_occ:
            for a, b in zip(plist, plist[1:]):
                counts[b - a] += 1
    return counts


def estimate_period(seq, min_period=5000, minimap2="minimap2", threads=1,
                    use_minimap2=True):
    """Return (period or None, method)."""
    L = len(seq)
    if L < 2 * min_period:
        return None, "too_short"
    shifts = minimap2_shifts(seq, minimap2, threads) if use_minimap2 else None
    if shifts is not None:
        good = [(s, b) for s, b in shifts if s >= min_period and b >= 0.5 * min(s, L - s)]
        if not good:
            return None, "minimap2"
        return min(s for s, _ in good), "minimap2"
    counts = kmer_shifts(seq)
    # bin distances (copies may differ by small indels) and require support
    binned = collections.Counter()
    for d, n in counts.items():
        if d >= min_period:
            binned[int(round(d / 50.0))] += n
    if not binned:
        return None, "kmer"
    b, support = binned.most_common(1)[0]
    if support < 0.2 * min_period:
        return None, "kmer"
    # exact mode within the winning bin
    cands = [(n, d) for d, n in counts.items() if int(round(d / 50.0)) == b]
    return max(cands)[1], "kmer"


def find_boundaries(seq, period, circular, k=K, tol=0.1):
    """Homologous cut points (one per copy) using an anchor k-mer."""
    L = len(seq)
    ext = seq + seq[:k - 1] if circular else seq
    index = collections.defaultdict(list)
    for i in range(len(ext) - k + 1):
        index[ext[i:i + k]].append(i)
    n_expected = int(round(L / float(period))) if circular else L // period
    best = None
    step = max(1, period // 200)
    for x in range(0, min(period, L - k), step):
        occ = index.get(ext[x:x + k], [])
        if len(occ) < 2:
            continue
        gaps = [b - a for a, b in zip(occ, occ[1:])]
        if circular:
            gaps.append(occ[0] + L - occ[-1])
        if all(abs(g - period) <= tol * period for g in gaps):
            if len(occ) == n_expected:
                return occ
            if best is None or len(occ) > len(best):
                best = occ
    return best


def split_copies(seq, circular=True, min_period=5000, minimap2="minimap2",
                 threads=1, use_minimap2=True):
    """Return (copies, info).  ``copies`` is a list of (start, end, sequence)."""
    L = len(seq)
    period, method = estimate_period(seq, min_period, minimap2, threads, use_minimap2)
    info = {"length": L, "circular": circular, "period": period, "method": method}
    if not period or (circular and L < 1.5 * period) or (not circular and L < 2 * period):
        info["n_copies"] = 1
        return [(0, L, seq)], info
    bounds = find_boundaries(seq, period, circular)
    if not bounds:
        # fall back to evenly spaced cuts
        n = int(round(L / float(period))) if circular else L // period
        bounds = [i * period for i in range(n)]
        info["boundaries"] = "even"
    else:
        info["boundaries"] = "anchor"
    copies = []
    for i, b in enumerate(bounds):
        if i + 1 < len(bounds):
            e = bounds[i + 1]
            copies.append((b, e, seq[b:e]))
        elif circular:
            e = bounds[0] + L
            copies.append((b, e, (seq + seq)[b:e]))
        elif L - b >= period * (1 - 0.02):
            # trailing copy of a linear contig, keep only if complete
            copies.append((b, b + period, seq[b:b + period]))
    info["n_copies"] = len(copies)
    info["copy_lengths"] = [e - s for s, e, _ in copies]
    return copies, info
