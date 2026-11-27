"""Minimal FASTA/FASTQ helpers (no Biopython needed)."""

import gzip
import io
import os

_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                      "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def open_text(path, mode="rt"):
    """Open a plain or gzip-compressed file (detected by magic bytes)."""
    if "r" in mode:
        with open(path, "rb") as fh:
            magic = fh.read(2)
        if magic == b"\x1f\x8b":
            return io.TextIOWrapper(gzip.open(path, "rb"))
        return open(path, "r")
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "wb"))
    return open(path, mode.replace("t", ""))


def read_fasta(path):
    """Return a list of (name, description, sequence) tuples.

    ``name`` is the first word of the header; ``description`` the remainder.
    """
    records = []
    name = desc = None
    chunks = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, desc, "".join(chunks)))
                header = line[1:].strip().split(None, 1)
                name = header[0] if header else ""
                desc = header[1] if len(header) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, desc, "".join(chunks)))
    return records


def write_fasta(records, path, width=60):
    """Write (name, sequence) or (name, description, sequence) tuples."""
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w") as out:
        for rec in records:
            if len(rec) == 3:
                name, desc, seq = rec
            else:
                (name, seq), desc = rec, ""
            out.write(">%s%s\n" % (name, (" " + desc) if desc else ""))
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def rotate(seq, offset):
    """Rotate a circular sequence so that position ``offset`` becomes 0."""
    if not seq:
        return seq
    offset %= len(seq)
    return seq[offset:] + seq[:offset]
