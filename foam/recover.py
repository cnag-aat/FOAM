"""Sequence-aware repair of a MITOS annotation.

MITOS (with ``--best``-like behaviour of its final selection) regularly drops
genes that its own intermediate searches found, especially in compact or
divergent mitogenomes (e.g. stylommatophoran snails):

* rRNAs whose fuzzy ends overlap neighbouring tRNAs -> recovered from
  ``ignored.mitos`` and trimmed back to the neighbours;
* short/divergent protein-coding genes (nad4l, atp8, nad6) for which BLAST
  hits exist in ``blast/prot`` but no gene model was built -> rebuilt as the
  complete ORF around the best hit;
* aberrant tRNAs below MiTFi's thresholds -> rescued with a targeted
  ``cmsearch --max`` of the missing tRNA's covariance model;
* a protein whose 5' end is overlapped by a (recovered) tRNA -> start moved
  to the next in-frame start codon that clears the tRNA, if one is close.

Every change is noted on the feature and in the change log.
"""

import glob
import os
import re
import shutil
import subprocess
import tempfile

from . import annotation as ann
from . import seqio
from .gencode import GeneticCode

_BLAST_RE = re.compile(r"-(\d+)\.([A-Za-z0-9]+)\.blast$")


# ------------------------------------------------------------------ helpers
def _sub(seq, start, end):
    """Circular substring [start, end)."""
    L = len(seq)
    start %= L
    end = start + (end - start) % L if end - start < L else start + L
    return (seq + seq)[start:end] if end > L else seq[start:end]


def _codon(seq, pos, strand):
    """Codon whose first (5') base is at ``pos`` on ``strand``."""
    L = len(seq)
    if strand > 0:
        return _sub(seq, pos % L, pos % L + 3)
    return seqio.revcomp(_sub(seq, (pos - 2) % L, (pos - 2) % L + 3))


def _overlaps(f, feats, L, types=None):
    out = []
    for g in feats:
        if g is f or (types and g.type not in types):
            continue
        o = ann.overlap(f, g, L)
        if o > 0:
            out.append((g, o))
    return out


# ------------------------------------------------------------------ rRNA
def trim_to_neighbours(f, feats, L, min_keep=0.7):
    """Trim a recovered RNA so it overlaps no other feature.

    A neighbour covering one end pushes that end past it; a neighbour lying
    wholly inside is cut away on whichever side is nearer.  Repeats until
    nothing overlaps, so the result does not depend on the order of
    ``feats`` (trimming past one gene can expose the next, e.g.
    cox1 / trnV / rrnL).  Fails if less than ``min_keep`` remains.
    """
    orig = f.length
    for _ in range(len(feats) + 1):
        changed = False
        for g, o in _overlaps(f, feats, L):
            if g.type == "rep_origin":
                continue
            for shift in (-L, 0, L):
                gs, ge = g.start + shift, g.end + shift
                if ge <= f.start or gs >= f.end:
                    continue
                if gs <= f.start:
                    f.start = ge
                elif ge >= f.end:
                    f.end = gs
                elif gs - f.start <= f.end - ge:
                    f.start = ge
                else:
                    f.end = gs
                changed = True
                break
            if changed or f.end - f.start < min_keep * orig:
                break
        if not changed or f.end - f.start < min_keep * orig:
            break
    if f.end - f.start < min_keep * orig:
        return False
    ann.normalize(f, L)
    if f.length != orig:
        f.notes.append("ends trimmed to neighbouring genes (%d -> %d bp)" % (orig, f.length))
    return True


def recover_rrnas(feats, candidates, expected, L):
    log = []
    present = set(f.name for f in feats)
    for gene in [g for g in expected if g.startswith("rrn") and g not in present]:
        cands = sorted([c for c in candidates if c.name == gene], key=ann._candidate_rank)
        for c in cands:
            nf = c.copy_feature()
            nf.copy = None
            nf.source = "recovered:" + c.source
            if not trim_to_neighbours(nf, feats, L):
                continue
            nf.notes.append("recovered from %s (MITOS dropped it for overlaps); check ends" % c.source)
            feats.append(nf)
            log.append("recovered %s at %d-%d from %s, trimmed to neighbours" % (
                gene, nf.start + 1, nf.end, c.source))
            break
    return log


# ------------------------------------------------------------------ proteins
def parse_blast_hits(mitos_dir, gene, L, max_evalue=1e-5):
    """MITOS protein BLAST hits for ``gene`` in genome coordinates.

    Returns list of dicts {start, end, strand, sstart, bits, evalue} with
    ``start`` = 5' end position (0-based, of the first codon base on that
    strand) and ``end`` = 3' end position.
    """
    hits = []
    for path in glob.glob(os.path.join(mitos_dir, "blast", "prot", "*.%s.blast" % gene)):
        m = _BLAST_RE.search(os.path.basename(path))
        if not m or m.group(2) != gene:
            continue
        off = int(m.group(1))
        with open(path) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 12:
                    continue
                ev = float(c[10])
                if ev > max_evalue:
                    continue
                qs, qe = int(c[6]) - 1, int(c[7]) - 1
                strand = 1 if qe >= qs else -1
                hits.append({"start": (qs + off) % L, "end": (qe + off) % L, "strand": strand,
                             "sstart": int(c[8]), "bits": float(c[11]), "evalue": ev})
    hits.sort(key=lambda h: -h["bits"])
    return hits


def orf_from_hit(seq, hit, code, max_upstream_codons=30, max_codons=1000):
    """Complete ORF (start..stop) in the frame of a protein hit."""
    gc = GeneticCode(code)
    L, st = len(seq), hit["strand"]
    step = 3 * st
    p5 = hit["start"]
    # 3' end: walk from the hit's 5' end to the first stop codon
    p, n = p5, 0
    while n < max_codons:
        if _codon(seq, p, st) in gc.stops:
            break
        p += step
        n += 1
    else:
        return None
    stop_pos = p
    # 5' end: walk upstream (no stops) collecting start codons, at most
    # (hit's offset in the reference protein + max_upstream_codons) codons;
    # take the most upstream one.
    starts, p, n = [], p5, 0
    while n < (hit["sstart"] - 1) + max_upstream_codons:
        if _codon(seq, p, st) in gc.starts:
            starts.append(p)
        q = p - step
        if _codon(seq, q, st) in gc.stops:
            break
        p, n = q, n + 1
    if not starts:
        return None
    # prefer codons read as Met (ATG, and ATA in most mito codes) over
    # alternative initiators (ATT, ATC, TTG, GTG)
    met = [s for s in starts if gc.codon[_codon(seq, s, st)] == "M"]
    start_pos = (met or starts)[-1]
    # convert to half-open + coordinates
    if st > 0:
        s, e = start_pos, stop_pos + 3
    else:
        s, e = stop_pos - 2, start_pos + 1
    f = ann.normalize(ann.Feature("?", "gene", s, e, st, method="orf"), L)
    cds = _sub(seq, f.start, f.end)
    if st < 0:
        cds = seqio.revcomp(cds)
    return f, gc.translate(cds)


def recover_proteins(feats, mitos_dir, seq, expected, code, L, max_overlap=120):
    log = []
    present = set(f.name for f in feats)
    prot_names = set(ann.METAZOA_PROT)
    for gene in [g for g in expected if g in prot_names and g not in present]:
        for hit in parse_blast_hits(mitos_dir, gene, L)[:5]:
            res = orf_from_hit(seq, hit, code)
            if res is None:
                continue
            f, prot = res
            if len(prot) < 30:
                continue
            ovs = [(g, o) for g, o in _overlaps(f, feats, L) if g.type != "rep_origin"]
            if any(o > max_overlap for g, o in ovs) or \
                    any(g.type == "gene" and o > 0.5 * min(g.length, f.length) for g, o in ovs):
                continue
            f.name = gene
            f.source = "recovered:blast/prot"
            f.score = hit["bits"]
            f.notes.append("rebuilt as complete ORF around MITOS blastx hit (bits %g, e %g); "
                           "%d aa" % (hit["bits"], hit["evalue"], len(prot) - 1))
            big = ["%s %d bp" % (g.name, o) for g, o in ovs if o > 15]
            if big:
                f.notes.append("overlaps " + ", ".join(big) + "; check boundaries")
            feats.append(f)
            log.append("recovered %s at %d-%d (%s) as ORF of %d aa from blast hits%s" % (
                gene, f.start + 1, f.end, "+" if f.strand > 0 else "-", len(prot) - 1,
                ("; overlaps " + ", ".join(big)) if big else ""))
            break
    return log


# ------------------------------------------------------------------ tRNAs
def _cm_for(refdir, refseq, name):
    p = os.path.join(refdir, refseq, "ncRNA", name + ".cm")
    return p if os.path.exists(p) else None


def cmsearch_trna(seq, cm, cmsearch="cmsearch", max_evalue=0.01, threads=1):
    """Hits of one tRNA model on a circular sequence (genome coordinates)."""
    exe = shutil.which(cmsearch) if os.path.sep not in cmsearch else cmsearch
    if not exe:
        return None
    L = len(seq)
    tmpd = tempfile.mkdtemp(prefix="foam_cm_")
    try:
        fa = os.path.join(tmpd, "s.fa")
        tbl = os.path.join(tmpd, "hits.tbl")
        seqio.write_fasta([("s", seq + seq[:200])], fa)
        rc = subprocess.call([exe, "--max", "-E", str(max_evalue), "--cpu", str(threads),
                              "--tblout", tbl, cm, fa],
                             stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        if rc != 0:
            return []
        hits = []
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.split()
                a, b, strand = int(c[7]), int(c[8]), c[9]
                lo, hi = min(a, b) - 1, max(a, b)
                hits.append({"start": lo % L, "end": lo % L + (hi - lo),
                             "strand": 1 if strand == "+" else -1,
                             "bits": float(c[14]), "evalue": float(c[15])})
        return hits
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def rescue_trnas(feats, seq, expected, L, refdir, refseq, cmsearch="cmsearch",
                 max_evalue=0.2, max_rna_overlap=10, max_pcg_overlap=40, threads=1):
    """Genome-wide e-values: a tRNA found at e~0.002 in a 300 bp gap scores
    e~0.1 on a 15 kb genome, hence the lenient default.  Hits on existing
    tRNAs/rRNAs (other isotypes) or deep inside proteins are rejected."""
    log = []
    present = set(f.name for f in feats)
    for gene in [g for g in expected if g.startswith("trn") and g not in present]:
        cm = _cm_for(refdir, refseq, gene)
        if not cm:
            continue
        hits = cmsearch_trna(seq, cm, cmsearch, max_evalue, threads)
        if hits is None:
            log.append("cmsearch not found; could not rescue %s" % gene)
            return log
        for h in sorted(hits, key=lambda h: h["evalue"]):
            f = ann.normalize(ann.Feature(gene, "tRNA", h["start"], h["end"], h["strand"],
                                          score=h["evalue"], bitscore=h["bits"],
                                          method="cmsearch", source="rescued:cmsearch"), L)
            ovs = _overlaps(f, feats, L)
            if any(g.type in ("tRNA", "rRNA") and o > max_rna_overlap for g, o in ovs):
                continue      # this is another tRNA hit by a similar model
            if any(g.type == "gene" and o > max_pcg_overlap for g, o in ovs):
                continue
            f.anticodon = _anticodon_guess(seq, f)
            f.notes.append("rescued by cmsearch --max with %s model (e-value %g); "
                           "below MITOS thresholds, please check" % (gene, h["evalue"]))
            feats.append(f)
            log.append("rescued %s at %d-%d (%s) by cmsearch, e-value %g" % (
                gene, f.start + 1, f.end, "+" if f.strand > 0 else "-", h["evalue"]))
            break
    return log


# Standard metazoan mt anticodons for isotypes that have a single one.
STD_ANTICODON = {"trnA": "TGC", "trnR": "TCG", "trnN": "GTT", "trnD": "GTC",
                 "trnC": "GCA", "trnQ": "TTG", "trnE": "TTC", "trnG": "TCC",
                 "trnH": "GTG", "trnI": "GAT", "trnK": "TTT", "trnM": "CAT",
                 "trnF": "GAA", "trnP": "TGG", "trnT": "TGT", "trnW": "TCA",
                 "trnY": "GTA", "trnV": "TAC", "trnL1": "TAG", "trnL2": "TAA",
                 "trnS2": "TGA"}


def _anticodon_guess(seq, f):
    """Expected anticodon if it occurs (preceded by T, i.e. U33) in the
    middle of the tRNA; otherwise None."""
    ac = STD_ANTICODON.get(f.name)
    if not ac:
        return None
    t = _sub(seq, f.start, f.end)
    if f.strand < 0:
        t = seqio.revcomp(t)
    lo, hi = int(len(t) * 0.3), int(len(t) * 0.7)
    return ac if ("T" + ac) in t[lo:hi] else None


# ------------------------------------------------------------------ starts
def fix_5prime_overlaps(feats, seq, code, L, min_overlap=20, max_shift_codons=20):
    """Move a protein's start downstream past an overlapping tRNA if an in-frame
    start codon exists within ``max_shift_codons``."""
    gc = GeneticCode(code)
    log = []
    for p in [f for f in feats if f.type == "gene" and f.name in ann.METAZOA_PROT]:
        five = p.start if p.strand > 0 else p.end - 1
        for t in [f for f in feats if f.type == "tRNA"]:
            o = ann.overlap(p, t, L)
            if o < min_overlap:
                continue
            # tRNA must cover the protein's 5' end
            probe = ann.Feature("x", "x", five % L, five % L + 1, 1)
            if ann.overlap(probe, t, L) == 0:
                continue
            step = 3 * p.strand
            pos = five
            for k in range(1, max_shift_codons + 1):
                pos = five + k * step
                probe = ann.Feature("x", "x", pos % L, pos % L + 1, 1)
                if ann.overlap(probe, t, L) == 0 and _codon(seq, pos, p.strand) in gc.starts:
                    old = (p.start + 1, p.end)
                    if p.strand > 0:
                        p.start = pos
                    else:
                        p.end = pos + 1
                    ann.normalize(p, L)
                    p.notes.append("start moved %d codons downstream to %s to clear %s" % (
                        k, _codon(seq, pos, p.strand), t.name))
                    log.append("%s start moved %d codons (%d-%d -> %d-%d) to clear %s" % (
                        p.name, k, old[0], old[1], p.start + 1, p.end, t.name))
                    break
            break
    return log
