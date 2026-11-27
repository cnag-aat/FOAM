#!/usr/bin/env python3
"""Stand-in for MITOS' runmitos for tests.

Finds marker genes (sequences in markers.fa next to this script, header =
'name type') on both strands of the circular input and writes result.mitos /
ignored.mitos like MITOS 2.1.x (0-based inclusive coordinates, stop < start
for features crossing the origin).
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from foam import seqio  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("-i")
p.add_argument("-c")
p.add_argument("-o")
p.add_argument("-r")
p.add_argument("-R")
p.add_argument("--noplots", action="store_true")
a, _ = p.parse_known_args()

markers = os.environ.get("FAKE_MITOS_MARKERS", os.path.join(HERE, "markers.fa"))
name, _d, seq = seqio.read_fasta(a.i)[0]
L = len(seq)
dbl = seq + seq
lines = []
for mname, mtype, mseq in seqio.read_fasta(markers):
    for strand, probe in ((1, mseq), (-1, seqio.revcomp(mseq))):
        pos = dbl.find(probe)
        if 0 <= pos < L:
            start, stop = pos, (pos + len(probe) - 1) % L
            ac = "GAA" if mname == "trnF" else "-"
            lines.append("\t".join(map(str, [name, mtype, mname, "mitfi" if mtype != "gene" else "mitos",
                                             start, stop, strand, 1e-10, 50.0, ac, "-", ".", ".", "."])))
os.makedirs(a.o, exist_ok=True)
with open(os.path.join(a.o, "result.mitos"), "w") as fh:
    fh.write("\n".join(lines) + ("\n" if lines else ""))
open(os.path.join(a.o, "ignored.mitos"), "w").close()
