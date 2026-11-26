"""Select mitochondrial reads by k-mer abundance (coverage), not identity.

Input is a *candidate* read set from a permissive prefilter (any read with a
few hundred bp aligned to a bait or to a mitochondrial protein), so it holds
the true mito reads plus reads carrying NUMTs or other partial matches.

Within that set, k-mers from the mitogenome are seen hundreds to thousands
of times (mito copy number), while k-mers from the nuclear flanks of a NUMT
read occur only in the handful of reads from that locus, and sequencing-error
k-mers occur once or twice.  So:

1. Count canonical k-mers (a fixed hash-sampled 1/``mod`` subset, which is
   the same subset in every read) over the candidate reads.
2. Find the mito coverage peak ``P``: the abundance carrying the most k-mer
   occurrences on a log scale.  ``T = coverage_fraction * P`` separates mito
   k-mers from nuclear/error k-mers.
3. Score each read: fraction of its sampled k-mers with abundance >= T.  A
   mito read scores the fraction of its k-mers that are error-free; a NUMT
   read scores only its mito-like stretch, diluted by its nuclear flanks.
4. The typical mito read score ``S`` is the mode of the score distribution
   among reads that score at least half the best scores; keep reads with
   score >= ``min_score_fraction * S``.

Both cut-offs are relative to peaks found in the data, so they do not
depend on the read accuracy or on how divergent the bait was.
"""

import collections
import math
import random

from . import seqio

_COMP = str.maketrans("ACGT", "TGCA")


def iter_reads(path):
    """(id, sequence) from FASTA or FASTQ, plain or gzipped."""
    with seqio.open_text(path) as fh:
        first = fh.readline()
        if not first:
            return
        if first.startswith(">"):
            rid, chunks = first[1:].split(None, 1)[0], []
            for line in fh:
                if line.startswith(">"):
                    yield rid, "".join(chunks)
                    rid, chunks = line[1:].split(None, 1)[0], []
                else:
                    chunks.append(line.strip())
            yield rid, "".join(chunks)
        else:
            line = first
            while line:
                rid = line[1:].split(None, 1)[0]
                seq = fh.readline().strip()
                fh.readline()
                fh.readline()
                yield rid, seq
                line = fh.readline()


def sampled_kmers(seq, k, mod):
    """Canonical k-mers whose hash is 0 modulo ``mod`` (hash is stable
    within one Python process, which is all we need)."""
    s = seq.upper()
    r = s.translate(_COMP)[::-1]
    n = len(s)
    out = []
    for i in range(n - k + 1):
        f = s[i:i + k]
        c = r[n - k - i:n - i]
        if c < f:
            f = c
        if hash(f) % mod == 0 and "N" not in f:
            out.append(f)
    return out


def count(paths, k=21, mod=8, max_bases=300000000, seed=1):
    """Count sampled k-mers.  Above ``max_bases`` reads are subsampled."""
    total = 0
    for p in paths:
        for _rid, s in iter_reads(p):
            total += len(s)
    frac = min(1.0, float(max_bases) / total) if total else 1.0
    rnd = random.Random(seed)
    counts = collections.Counter()
    used = 0
    for p in paths:
        for _rid, s in iter_reads(p):
            if frac < 1.0 and rnd.random() > frac:
                continue
            counts.update(sampled_kmers(s, k, mod))
            used += len(s)
    # abundances from a subsample are scaled back to the full set
    return counts, {"total_bases": total, "counted_bases": used, "sample_fraction": frac}


def coverage_peak(counts, min_abundance=5, bins_per_doubling=4):
    """Abundance at the occurrence-weighted mode (log-binned), or None."""
    weight = collections.Counter()
    members = collections.defaultdict(list)
    for a in counts.values():
        if a < min_abundance:
            continue
        b = int(math.log(a, 2) * bins_per_doubling)
        weight[b] += a
        members[b].append(a)
    if not weight:
        return None
    b = max(weight, key=lambda x: (weight[x], x))
    vals = sorted(members[b])
    return vals[len(vals) // 2]


def histogram(counts, bins_per_doubling=2):
    h = collections.Counter()
    for a in counts.values():
        h[2 ** (int(math.log(a, 2) * bins_per_doubling) / float(bins_per_doubling))] += a
    return [[round(k, 1), v] for k, v in sorted(h.items())]


def typical_score(scores, bin_width=0.02):
    """Mode of read scores among the upper part of the distribution."""
    if not scores:
        return None
    srt = sorted(scores)
    top = srt[int(0.95 * (len(srt) - 1))]
    upper = [s for s in srt if s >= 0.5 * top]
    h = collections.Counter(int(s / bin_width) for s in upper)
    b = max(h, key=lambda x: (h[x], x))
    return (b + 0.5) * bin_width


def select(paths, k=21, mod=8, coverage_fraction=0.2, min_score_fraction=0.5,
           min_kmers=20, max_bases=300000000):
    """Return (kept read ids, report dict)."""
    counts, info = count(paths, k, mod, max_bases)
    peak = coverage_peak(counts)
    report = dict(info, k=k, sampling=mod, coverage_peak=peak,
                  kmer_histogram=histogram(counts))
    if peak is None:
        report["error"] = "no high-abundance k-mers: no mitochondrial reads in the candidate set?"
        return [], report
    # peak and counts are on the same (possibly subsampled) scale
    thr = max(3, coverage_fraction * peak)
    report["abundance_threshold"] = thr
    ids, scores = [], []
    too_short = 0
    for p in paths:
        for rid, s in iter_reads(p):
            km = sampled_kmers(s, k, mod)
            if len(km) < min_kmers:
                too_short += 1
                continue
            sc = sum(1 for x in km if counts.get(x, 0) >= thr) / float(len(km))
            ids.append(rid)
            scores.append(sc)
    S = typical_score(scores)
    report["typical_mito_read_score"] = S
    report["score_threshold"] = None if S is None else min_score_fraction * S
    keep = [] if S is None else [r for r, sc in zip(ids, scores) if sc >= min_score_fraction * S]
    report.update(n_candidates=len(ids) + too_short, n_too_short=too_short, n_kept=len(keep),
                  score_histogram=sorted(collections.Counter(round(x, 1) for x in scores).items()))
    return keep, report
