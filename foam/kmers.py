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

Counting uses FastK (Gene Myers) when it is on PATH: all k-mers, multi-
threaded, with per-read count profiles.  The pure-Python fallback samples
1/``mod`` of the k-mers and is ~100x slower.
"""

import collections
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from multiprocessing.pool import ThreadPool

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
    """Abundance at the occurrence-weighted mode (log-binned), or None.
    ``counts``: k-mer -> abundance."""
    return peak_from_hist(collections.Counter(counts.values()), min_abundance,
                          bins_per_doubling)


def peak_from_hist(hist, min_abundance=5, bins_per_doubling=4):
    """Same as coverage_peak, from a histogram {abundance: distinct k-mers}."""
    weight = collections.Counter()
    members = collections.defaultdict(collections.Counter)
    for a, n in hist.items():
        if a < min_abundance or n <= 0:
            continue
        b = int(math.log(a, 2) * bins_per_doubling)
        weight[b] += a * n
        members[b][a] += n
    if not weight:
        return None
    b = max(weight, key=lambda x: (weight[x], x))
    m = members[b]
    half, cum = sum(m.values()) // 2, 0
    for a in sorted(m):
        cum += m[a]
        if cum > half:
            return a


def histogram(counts, bins_per_doubling=2):
    return histogram_from_hist(collections.Counter(counts.values()), bins_per_doubling)


def histogram_from_hist(hist, bins_per_doubling=2):
    """Log-binned k-mer occurrences per abundance, for the report."""
    h = collections.Counter()
    for a, n in hist.items():
        if a > 0 and n > 0:
            h[2 ** (int(math.log(a, 2) * bins_per_doubling) / float(bins_per_doubling))] += a * n
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


# --------------------------------------------------------------- FastK backend
FASTK_TOOLS = ("FastK", "Histex", "Profex")


def fastk_available(bindir=None):
    return all(shutil.which(os.path.join(bindir, t) if bindir else t) for t in FASTK_TOOLS)


class FastK(object):
    """Count all canonical k-mers with Gene Myers' FastK and stream each
    read's k-mer count profile (``Profex -zA``: runs of equal count)."""

    def __init__(self, paths, k=21, threads=4, memory=4, tmpdir=None, bindir=None):
        self.bin = (lambda t: os.path.join(bindir, t)) if bindir else (lambda t: t)
        self.dir = tempfile.mkdtemp(prefix="foam_fastk.", dir=tmpdir)
        self.root = os.path.join(self.dir, "kmers")
        # FastK drops some records (e.g. empty FASTA entries) and profiles are
        # identified only by ordinal, so feed it one FASTA of the reads with
        # >= k bases and keep their names here, in the same order.
        self.ids, self.short_ids = [], []
        fa = os.path.join(self.dir, "reads.fa")
        with open(fa, "w") as out:
            for p in paths:
                for rid, s in iter_reads(p):
                    if len(s) >= k:
                        self.ids.append(rid)
                        out.write(">%d\n%s\n" % (len(self.ids), s))
                    else:
                        self.short_ids.append(rid)
        # FastK 1.2 sometimes fails to open its own partition files at some
        # thread counts (seen with -T6..8); retry with fewer threads.
        errors = []
        for t in sorted({max(1, threads), min(threads, 4), 1}, reverse=True):
            cmd = [self.bin("FastK"), "-k%d" % k, "-p", "-T%d" % t, "-M%d" % memory,
                   "-N" + self.root, "-P" + self.dir, fa]
            r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
            if r.returncode == 0:
                self.threads = t
                break
            errors.append("%s\n%s" % (" ".join(cmd), r.stderr.strip()[:2000]))
            for f in os.listdir(self.dir):
                if f != "reads.fa":
                    os.remove(os.path.join(self.dir, f))
        os.remove(fa)
        if r.returncode != 0:
            self.close()
            raise RuntimeError("FastK failed:\n" + "\n".join(errors))

    def hist(self):
        """{abundance: distinct k-mers}"""
        out = subprocess.check_output([self.bin("Histex"), "-A", "-h1:32767", self.root],
                                      universal_newlines=True)
        h = {}
        for line in out.splitlines():
            f = line.split()
            if len(f) == 2 and f[0].isdigit() and f[1].isdigit():
                h[int(f[0])] = int(f[1])
        return h

    def scores(self, thr, workers=1):
        """Yield (read id, n_kmers, n_kmers with count >= thr) in input order.
        Profex is single-threaded, so read ranges go through parallel
        Profex | reducer pipelines (threads here only wait on them)."""
        n = len(self.ids)
        if n == 0:
            return
        w = max(1, min(workers, n // 500 or 1))
        step = -(-n // w)
        jobs = [(self.bin("Profex"), self.root, a + 1, min(n, a + step), thr)
                for a in range(0, n, step)]
        if w == 1:
            results = map(_profex_range, jobs)
        else:
            pool = ThreadPool(w)
            results = pool.imap(_profex_range, jobs)
        i = 0
        try:
            for (first, last), rows in results:
                if len(rows) != last - first + 1:
                    raise RuntimeError("Profex returned %d profiles for reads %d-%d"
                                       % (len(rows), first, last))
                for tot, hi in rows:
                    yield self.ids[i], tot, hi
                    i += 1
        finally:
            if w > 1:
                pool.terminate()

    def close(self):
        shutil.rmtree(self.dir, ignore_errors=True)


# Reduces ``Profex -zA`` output (Read <n> / start end count lines) to one
# "n_kmers<TAB>n_kmers_at_or_above_thr" line per read.  Run as its own
# process so that several ranges are parsed on several cores.
_REDUCER = r"""
import sys
thr = float(sys.argv[1])
out, tot, hi = [], -1, 0
for line in sys.stdin.buffer:
    f = line.split(b"\t")
    if len(f) == 3:
        ln = int(f[1]) - int(f[0]) + 1
        tot += ln
        if int(f[2]) >= thr:
            hi += ln
    elif f[0] == b"Read":
        if tot >= 0:
            out.append("%d\t%d" % (tot, hi))
        tot, hi = 0, 0
if tot >= 0:
    out.append("%d\t%d" % (tot, hi))
sys.stdout.write("\n".join(out) + ("\n" if out else ""))
"""


def _profex_range(job):
    """[(n_kmers, n_kmers >= thr)] for reads first..last (1-based)."""
    profex, root, first, last, thr = job
    p1 = subprocess.Popen([profex, "-zA", root, "%d-%d" % (first, last)],
                          stdout=subprocess.PIPE)
    p2 = subprocess.Popen([sys.executable, "-c", _REDUCER, repr(float(thr))],
                          stdin=p1.stdout, stdout=subprocess.PIPE, universal_newlines=True)
    p1.stdout.close()
    out = p2.communicate()[0]
    if p1.wait() != 0 or p2.returncode != 0:
        raise RuntimeError("Profex failed on %s (reads %d-%d)" % (root, first, last))
    rows = [tuple(int(x) for x in line.split("\t")) for line in out.splitlines()]
    return (first, last), rows


# ------------------------------------------------------------------ selection
def select(paths, k=21, mod=8, coverage_fraction=0.2, min_score_fraction=0.5,
           min_kmers=20, max_bases=300000000, backend="auto", threads=4, memory=4,
           tmpdir=None, fastk_bin=None):
    """Return (kept read ids, report dict).

    backend: ``fastk`` (all k-mers, FastK; fast, multi-threaded), ``python``
    (hash-sampled 1/``mod`` of the k-mers, subsampled above ``max_bases``;
    no dependencies) or ``auto`` (FastK when on PATH).  ``min_kmers`` is in
    sampled k-mers, i.e. a read needs about ``min_kmers * mod`` k-mers.
    """
    if backend == "auto":
        backend = "fastk" if fastk_available(fastk_bin) else "python"
    if backend == "fastk":
        return _select_fastk(paths, k, coverage_fraction, min_score_fraction,
                             min_kmers * mod, threads, memory, tmpdir, fastk_bin)
    if backend != "python":
        raise ValueError("unknown k-mer backend %r" % backend)
    counts, info = count(paths, k, mod, max_bases)
    peak = coverage_peak(counts)
    report = dict(info, backend="python", k=k, sampling=mod, coverage_peak=peak,
                  kmer_histogram=histogram(counts))
    if peak is None:
        return [], _no_peak(report)
    # peak and counts are on the same (possibly subsampled) scale
    thr = max(3, coverage_fraction * peak)
    report["abundance_threshold"] = thr

    def scored():
        for p in paths:
            for rid, s in iter_reads(p):
                km = sampled_kmers(s, k, mod)
                if len(km) < min_kmers:
                    yield rid, None
                else:
                    yield rid, sum(1 for x in km if counts.get(x, 0) >= thr) / float(len(km))
    return _finish(scored(), min_score_fraction, report)


def _select_fastk(paths, k, coverage_fraction, min_score_fraction, min_kmers,
                  threads, memory, tmpdir, fastk_bin):
    t0 = time.time()
    fk = FastK(paths, k, threads, memory, tmpdir, fastk_bin)
    try:
        hist = fk.hist()
        peak = peak_from_hist(hist)
        report = dict(backend="fastk", fastk_threads=fk.threads, k=k, sampling=1,
                      coverage_peak=peak,
                      total_kmers=sum(a * n for a, n in hist.items()),
                      kmer_histogram=histogram_from_hist(hist))
        if peak is None:
            return [], _no_peak(report)
        thr = max(3, coverage_fraction * peak)
        report["abundance_threshold"] = thr

        def scored():
            for rid in fk.short_ids:
                yield rid, None
            for rid, n, hi in fk.scores(thr, threads):
                yield rid, (None if n < min_kmers else hi / float(n))
        keep, report = _finish(scored(), min_score_fraction, report)
        report["seconds"] = round(time.time() - t0, 1)
        return keep, report
    finally:
        fk.close()


def _no_peak(report):
    report["error"] = "no high-abundance k-mers: no mitochondrial reads in the candidate set?"
    return report


def _finish(scored, min_score_fraction, report):
    ids, scores, too_short = [], [], 0
    for rid, sc in scored:
        if sc is None:
            too_short += 1
        else:
            ids.append(rid)
            scores.append(sc)
    S = typical_score(scores)
    report["typical_mito_read_score"] = S
    report["score_threshold"] = None if S is None else min_score_fraction * S
    keep = [] if S is None else [r for r, sc in zip(ids, scores) if sc >= min_score_fraction * S]
    report.update(n_candidates=len(ids) + too_short, n_too_short=too_short, n_kept=len(keep),
                  score_histogram=sorted(collections.Counter(round(x, 1) for x in scores).items()))
    return keep, report
