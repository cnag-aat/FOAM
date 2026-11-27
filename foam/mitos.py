"""Run MITOS2 and turn its output into a repaired, scored annotation."""

import logging
import os
import shlex
import shutil
import subprocess

from . import annotation as ann
from . import recover, seqio

log = logging.getLogger("foam")


class MitosError(RuntimeError):
    pass


_MITOS_NAMES = ("runmitos", "runmitos.py")  # >= 2.1.7, <= 2.1.6


def mitos_command(runmitos="auto"):
    """The MITOS command as a list.

    ``auto`` (or a bare ``runmitos``/``runmitos.py`` that is not on PATH)
    resolves to whichever of the two names is installed, so one config works
    with MITOS 2.1.3 and 2.1.10.
    """
    words = shlex.split(runmitos or "auto")
    if words == ["auto"] or (len(words) == 1 and words[0] in _MITOS_NAMES
                             and not shutil.which(words[0])):
        for name in _MITOS_NAMES:
            if shutil.which(name):
                return [name]
        raise MitosError("neither runmitos nor runmitos.py is on PATH")
    return words


def run_mitos(fasta, outdir, code, refdir, refseq, runmitos="auto",
              extra_opts="", logfile=None):
    """Run MITOS (circular mode, no plots) into ``outdir``."""
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    cmd = mitos_command(runmitos) + [
        "-i", fasta, "-c", str(code), "-o", outdir,
        "-r", refseq, "-R", refdir, "--noplots"] + shlex.split(extra_opts or "")
    log.info("running: %s", " ".join(shlex.quote(c) for c in cmd))
    with open(logfile or os.path.join(outdir, "mitos.log"), "w") as lf:
        rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0:
        raise MitosError("MITOS failed (exit %d) on %s; see %s" % (
            rc, fasta, logfile or os.path.join(outdir, "mitos.log")))
    return outdir


def annotate_dir(mitos_dir, seq, expected, code, refdir=None, refseq=None,
                 cmsearch="cmsearch", max_overlap=15, max_evalue=0.01,
                 keep_origins=False, rescue=True, threads=1):
    """Parse a MITOS run, repair it and return (features, qc, changelog)."""
    length = len(seq)
    feats = ann.parse_mitos_file(os.path.join(mitos_dir, "result.mitos"), length)
    if not keep_origins:
        # MITOS' OH/OL calls are BLAST/CM guesses and unreliable outside vertebrates
        feats = [f for f in feats if f.type != "rep_origin"]
    changes = ann.relabel_leu_ser(feats)
    cands = ann.parse_mitos_file(os.path.join(mitos_dir, "ignored.mitos"), length)
    for c in cands:
        c.source = "ignored.mitos"
    cands += ann.parse_mitfi_dir(mitos_dir, length)
    if rescue:
        changes += recover.recover_rrnas(feats, cands, expected, length)
    changes += ann.recover_missing(feats, cands, expected, length,
                                   max_overlap=max_overlap, max_evalue=max_evalue)
    if rescue:
        changes += recover.recover_proteins(feats, mitos_dir, seq, expected, code, length)
        if refdir and refseq:
            changes += recover.rescue_trnas(feats, seq, expected, length, refdir, refseq,
                                            cmsearch, threads=threads)
        changes += recover.fix_5prime_overlaps(feats, seq, code, length)
    feats.sort(key=lambda f: f.start)
    return feats, ann.qc(feats, expected, length), changes


def orient(seq, feats, anchors):
    """Rotate/reverse-complement so the first available anchor starts at 1 on +."""
    a = ann.choose_anchor(feats, anchors)
    if a is None:
        return seq, feats, None
    reverse, offset = ann.orientation_for(a, len(seq))
    s = seqio.revcomp(seq) if reverse else seq
    s = seqio.rotate(s, offset)
    return s, ann.transform(feats, len(seq), offset, reverse), {
        "anchor": a.name, "reverse_complemented": reverse, "rotation": offset}
