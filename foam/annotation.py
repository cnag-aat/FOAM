"""Mitogenome annotation model: parse MITOS output, QC, repair, orient, write.

Coordinates
-----------
A :class:`Feature` uses 0-based, half-open coordinates on a *circular*
sequence of length ``L``: ``0 <= start < L`` and ``start < end <= start + L``.
A feature with ``end > L`` wraps around the origin.

MITOS 2.1.x output used here
----------------------------
``result.mitos``   final annotation, 14 tab-separated columns
                   (acc, type, name, method, start, stop, strand, score,
                   bitscore, anticodon, anticodonpos, part, copy, struct);
                   start/stop are 0-based *inclusive*, strand is 1/-1.
``ignored.mitos``  same format; candidates MITOS dropped (overlaps etc.).
``mitfi-global/`` and ``mitfi-local/``  raw Infernal/MiTFi hits, files named
                   ``<fasta>-<offset>_<tRNA|rRNA|rep_origin>out.nc``; three
                   lines per hit, header fields separated by ``|``.  MITOS runs
                   on two rotations of a circular input, so the coordinates
                   are relative to the sequence rotated by ``<offset>``.
"""

import glob
import json
import math
import os
import re

# --------------------------------------------------------------------------
# Gene sets
# --------------------------------------------------------------------------

METAZOA_PROT = ["cox1", "cox2", "cox3", "cob", "atp6", "atp8", "nad1", "nad2",
                "nad3", "nad4", "nad4l", "nad5", "nad6"]
METAZOA_RRNA = ["rrnS", "rrnL"]
METAZOA_TRNA = ["trnA", "trnR", "trnN", "trnD", "trnC", "trnQ", "trnE", "trnG",
                "trnH", "trnI", "trnL1", "trnL2", "trnK", "trnM", "trnF",
                "trnP", "trnS1", "trnS2", "trnT", "trnW", "trnY", "trnV"]
GENE_SETS = {
    "metazoa": METAZOA_PROT + METAZOA_RRNA + METAZOA_TRNA,
    # Many cnidarians/sponges keep only trnM and trnW.
    "cnidaria": METAZOA_PROT + METAZOA_RRNA + ["trnM", "trnW"],
}

PRODUCTS = {
    "cox1": "cytochrome c oxidase subunit I",
    "cox2": "cytochrome c oxidase subunit II",
    "cox3": "cytochrome c oxidase subunit III",
    "cob": "cytochrome b",
    "atp6": "ATP synthase F0 subunit 6",
    "atp8": "ATP synthase F0 subunit 8",
    "nad1": "NADH dehydrogenase subunit 1",
    "nad2": "NADH dehydrogenase subunit 2",
    "nad3": "NADH dehydrogenase subunit 3",
    "nad4": "NADH dehydrogenase subunit 4",
    "nad4l": "NADH dehydrogenase subunit 4L",
    "nad5": "NADH dehydrogenase subunit 5",
    "nad6": "NADH dehydrogenase subunit 6",
    "rrnS": "12S ribosomal RNA",
    "rrnL": "16S ribosomal RNA",
}
TRNA_AA = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
           "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
           "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
           "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}


# Taxon presets: mitochondrial genetic code (NCBI transl_table), the gene whose
# 5' end becomes position 1, and the expected gene set.  trnF is the
# convention for vertebrates only; all other taxa start at the 5' end of cox1.
TAXON_PRESETS = {
    "vertebrata":    {"genetic_code": 2,  "anchor": ["trnF", "cox1"], "genes": "metazoa"},
    "mollusca":      {"genetic_code": 5,  "anchor": ["cox1"], "genes": "metazoa"},
    "arthropoda":    {"genetic_code": 5,  "anchor": ["cox1"], "genes": "metazoa"},
    "annelida":      {"genetic_code": 5,  "anchor": ["cox1"], "genes": "metazoa"},
    "nematoda":      {"genetic_code": 5,  "anchor": ["cox1"], "genes": "metazoa"},
    "invertebrate":  {"genetic_code": 5,  "anchor": ["cox1"], "genes": "metazoa"},
    "echinodermata": {"genetic_code": 9,  "anchor": ["cox1"], "genes": "metazoa"},
    "platyhelminthes": {"genetic_code": 9, "anchor": ["cox1"], "genes": "metazoa"},
    "tunicata":      {"genetic_code": 13, "anchor": ["cox1"], "genes": "metazoa"},
    "cnidaria":      {"genetic_code": 4,  "anchor": ["cox1"], "genes": "cnidaria"},
    "porifera":      {"genetic_code": 4,  "anchor": ["cox1"], "genes": "cnidaria"},
}
DEFAULT_ANCHOR = ["cox1"]


def taxon_preset(name):
    key = (name or "").strip().lower()
    if key not in TAXON_PRESETS:
        raise ValueError("unknown taxon group %r; choose from: %s" % (
            name, ", ".join(sorted(TAXON_PRESETS))))
    return TAXON_PRESETS[key]


def expected_genes(spec):
    """``spec`` is a gene-set name, a comma-separated string or a list."""
    if spec is None:
        spec = "metazoa"
    if isinstance(spec, (list, tuple)):
        return list(spec)
    if spec in GENE_SETS:
        return list(GENE_SETS[spec])
    return [g.strip() for g in spec.split(",") if g.strip()]


# --------------------------------------------------------------------------
# Feature
# --------------------------------------------------------------------------

class Feature(object):
    __slots__ = ("name", "type", "start", "end", "strand", "score",
                 "bitscore", "anticodon", "part", "copy", "method", "source",
                 "notes")

    def __init__(self, name, type, start, end, strand, score=None,
                 bitscore=None, anticodon=None, part=None, copy=None,
                 method="mitos", source="mitos", notes=None):
        self.name = name
        self.type = type
        self.start = start
        self.end = end
        self.strand = strand
        self.score = score
        self.bitscore = bitscore
        self.anticodon = anticodon
        self.part = part
        self.copy = copy
        self.method = method
        self.source = source
        self.notes = list(notes or [])

    def __repr__(self):
        return "Feature(%s %s %d-%d %s)" % (self.name, self.type, self.start,
                                            self.end, "+" if self.strand > 0 else "-")

    @property
    def length(self):
        return self.end - self.start

    def copy_feature(self):
        f = Feature(self.name, self.type, self.start, self.end, self.strand)
        for s in self.__slots__:
            setattr(f, s, getattr(self, s))
        f.notes = list(self.notes)
        return f

    def label(self):
        """MITOS-style display name, e.g. ``nad5_1-b`` or ``trnF(gaa)``."""
        s = self.name
        if self.copy is not None:
            s += "_%d" % self.copy
        if self.part is not None:
            s += "-" + chr(97 + self.part)
        if self.type == "tRNA":
            s += "(%s)" % (self.anticodon.lower() if self.anticodon else "---")
        return s


def normalize(f, length):
    """Put a feature into canonical circular coordinates (in place)."""
    span = f.end - f.start
    if span <= 0:
        span += length
    f.start %= length
    f.end = f.start + span
    return f


def overlap(a, b, length):
    """Overlap (nt) between two features on a circular sequence."""
    best = 0
    for shift in (-length, 0, length):
        lo = max(a.start, b.start + shift)
        hi = min(a.end, b.end + shift)
        best = max(best, hi - lo)
    return best


def _none(x):
    return None if x in (None, "", ".", "-", "None") else x


def _feature_from_mitos_cols(cols, length):
    if len(cols) == 14:
        acc, tpe, name, method, start, stop, strand, score, bits, ac, acpos, part, copy, struct = cols
    elif len(cols) == 13:
        acc, tpe, name, method, start, stop, strand, score, ac, part, copy, struct, acpos = cols
        bits = None
    else:
        raise ValueError("unexpected .mitos line with %d columns" % len(cols))
    start, stop, strand = int(start), int(stop), int(strand)
    if tpe == "tRNA" and not name.startswith("trn"):
        name = "trn" + name
    end = stop + 1
    if stop < start:        # wraps the origin
        end += length
    f = Feature(name=name, type=tpe, start=start, end=end,
                strand=1 if strand > 0 else -1,
                score=float(score) if _none(score) else None,
                bitscore=float(bits) if _none(bits) else None,
                anticodon=_none(ac),
                part=int(part) if _none(part) else None,
                copy=int(copy) if _none(copy) else None,
                method=method)
    return normalize(f, length)


def parse_mitos_file(path, length):
    feats = []
    if not os.path.exists(path):
        return feats
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            feats.append(_feature_from_mitos_cols(line.split("\t") if "\t" in line
                                                  else line.split(), length))
    return feats


_NC_RE = re.compile(r"-(\d+)_(tRNA|rRNA|rep_origin|intron)out\.nc$")


def parse_mitfi_dir(mitos_dir, length):
    """Raw MiTFi (Infernal) tRNA/rRNA hits from a MITOS output directory.

    Returns Features with ``score`` = e-value.  Coordinates are corrected for
    the rotation MITOS applied to the input.
    """
    hits = []
    for sub in ("mitfi-global", "mitfi-local"):
        for path in sorted(glob.glob(os.path.join(mitos_dir, sub, "*out.nc"))):
            m = _NC_RE.search(os.path.basename(path))
            offset = int(m.group(1)) if m else 0
            with open(path) as fh:
                lines = [l.strip() for l in fh if l.strip()]
            for i in range(0, len(lines) - 2, 3):
                head = lines[i]
                if not head.startswith(">"):
                    continue
                cols = head[1:].split("|")
                if len(cols) < 18:
                    continue
                start = int(cols[2]) - 1
                stop = int(cols[3]) - 1
                strand = 1 if int(cols[4]) > 0 else -1
                bits = float(cols[7])
                evalue = float(cols[8])
                anticodon = cols[12] if int(cols[11]) > 0 else None
                model = cols[14]
                if model.startswith("rrnS"):
                    tpe, name = "rRNA", "rrnS"
                elif model.startswith("rrnL"):
                    tpe, name = "rRNA", "rrnL"
                elif model.startswith("OL"):
                    tpe, name = "rep_origin", "OL"
                elif model.startswith("Intron"):
                    continue
                else:
                    mt = re.match(r"trn(\w\d*)", cols[13])
                    tpe = "tRNA"
                    name = "trn" + (mt.group(1) if mt else cols[13])
                lo, hi = min(start, stop), max(start, stop)
                f = Feature(name=name, type=tpe, start=lo + offset,
                            end=hi + 1 + offset, strand=strand, score=evalue,
                            bitscore=bits, anticodon=anticodon,
                            method="mitfi", source="raw:%s/%s" % (sub, os.path.basename(path)))
                hits.append(normalize(f, length))
    return hits


# --------------------------------------------------------------------------
# Repair / QC
# --------------------------------------------------------------------------

# Anticodon (5'->3') -> isotype for the two leucine and serine tRNAs.
# trnL1 decodes CUN (anticodon NAG), trnL2 decodes UUR (anticodon YAA).
# trnS1 decodes AGN (anticodon NCU), trnS2 decodes UCN (anticodon NGA).
def isotype_from_anticodon(name, anticodon):
    if not anticodon:
        return name
    ac = anticodon.upper().replace("U", "T")
    if len(ac) != 3:
        return name
    if name in ("trnL", "trnL1", "trnL2"):
        if ac[1:] == "AG":
            return "trnL1"
        if ac[1:] == "AA":
            return "trnL2"
    if name in ("trnS", "trnS1", "trnS2"):
        if ac[1:] == "CT":
            return "trnS1"
        if ac[1:] == "GA":
            return "trnS2"
    return name


def relabel_leu_ser(features):
    """Fix trnL/trnS isotype names using the anticodon.  Returns change log."""
    log = []
    for f in features:
        if f.type != "tRNA":
            continue
        new = isotype_from_anticodon(f.name, f.anticodon)
        if new != f.name:
            log.append("renamed %s(%s) at %d to %s" % (f.name, f.anticodon, f.start + 1, new))
            f.notes.append("renamed from %s by anticodon" % f.name)
            f.name = new
    return log


def _candidate_rank(f):
    """Smaller is better."""
    if f.type in ("tRNA", "rRNA", "rep_origin") and f.score is not None and f.method != "blast":
        return (0, f.score)                          # e-value
    return (1, -(f.score or 0.0))                    # MITOS protein quality


def recover_missing(features, candidates, expected, length,
                    max_overlap=15, max_evalue=0.01):
    """Add the best non-conflicting candidate for each missing expected gene.

    ``candidates`` are features from ``ignored.mitos`` and raw MiTFi hits.
    Returns a change log.
    """
    present = set(f.name for f in features)
    log = []
    for gene in expected:
        if gene in present:
            continue
        cands = [c for c in candidates if c.name == gene or
                 (c.type == "tRNA" and isotype_from_anticodon(c.name, c.anticodon) == gene)]
        if gene.startswith("trn") or gene.startswith("rrn"):
            cands = [c for c in cands if c.score is None or c.method != "mitfi" or c.score <= max_evalue]
        for c in sorted(cands, key=_candidate_rank):
            worst = max([overlap(c, f, length) for f in features
                         if f.type != "rep_origin"] or [0])
            if worst > max_overlap:
                continue
            nf = c.copy_feature()
            nf.name = gene
            nf.copy = None
            nf.source = "recovered:" + c.source
            ev = (" e-value %g" % c.score) if c.method == "mitfi" and c.score is not None else ""
            nf.notes.append("recovered from %s%s; please check" % (c.source, ev))
            features.append(nf)
            present.add(gene)
            log.append("recovered %s at %d-%d (%s)%s" % (gene, nf.start + 1, nf.end,
                                                         c.source, ev))
            break
    features.sort(key=lambda f: f.start)
    return log


def qc(features, expected, length):
    names = [f.name for f in features]
    present = [g for g in expected if g in names]
    missing = [g for g in expected if g not in names]
    fragmented = sorted(set(f.name for f in features if f.part is not None
                            and f.name in expected))
    duplicated = sorted(set(g for g in expected if
                            len(set((f.start, f.copy) for f in features
                                    if f.name == g and f.part in (None, 0))) > 1))
    prot_len = sum(f.length for f in features if f.type == "gene" and f.name in METAZOA_PROT)
    recovered = [f.name for f in features if f.source.startswith("recovered")]
    return {
        "length": length,
        "n_expected": len(expected),
        "n_present": len(present),
        "missing": missing,
        "fragmented": fragmented,
        "duplicated": duplicated,
        "recovered": recovered,
        "protein_coding_bp": prot_len,
    }


def completeness_key(q):
    """Sort key (bigger is better) for choosing among candidate copies."""
    return (q["n_present"] - len(q["recovered"]) * 0.5,
            -len(q["fragmented"]),
            -len(q["duplicated"]),
            q["protein_coding_bp"])


# --------------------------------------------------------------------------
# Orientation
# --------------------------------------------------------------------------

def transform(features, length, offset, reverse):
    """Coordinates after ``revcomp`` (if ``reverse``) then ``rotate(offset)``."""
    out = []
    for f in features:
        g = f.copy_feature()
        span = g.length
        if reverse:
            g.start = (length - f.end) % length
            g.strand = -f.strand
        g.start = (g.start - offset) % length
        g.end = g.start + span
        out.append(g)
    out.sort(key=lambda x: x.start)
    return out


def choose_anchor(features, anchors=("cox1",)):
    """Best copy of the first anchor gene present.  For a gene MITOS split
    into parts, part 0 ('-a') carries the 5' end."""
    for a in anchors:
        hits = [f for f in features if f.name == a and f.part in (None, 0)]
        if hits:
            hits.sort(key=_candidate_rank)
            return hits[0]
    return None


def orientation_for(anchor, length):
    """(reverse, offset) that puts the 5' end of ``anchor`` at position 0 on
    the + strand (for a - strand gene the 5' end is its highest coordinate)."""
    if anchor.strand > 0:
        return False, anchor.start
    return True, (length - anchor.end) % length


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def _gff_attr(d):
    def esc(v):
        return str(v).replace(";", "%3B").replace("=", "%3D").replace(",", "%2C")
    return ";".join("%s=%s" % (k, esc(v)) for k, v in d if v is not None and v != "")


def write_gff3(features, seqid, length, path, prefix, transl_table):
    lines = ["##gff-version 3", "##sequence-region %s 1 %d" % (seqid, length),
             "%s\tFOAM\tregion\t1\t%d\t.\t+\t.\t%s" % (
                 seqid, length, _gff_attr([("ID", seqid), ("Is_circular", "true"),
                                           ("genome", "mitochondrion"),
                                           ("transl_table", transl_table)]))]
    n = 0
    for f in sorted(features, key=lambda x: x.start):
        strand = "+" if f.strand > 0 else "-"
        s, e = f.start + 1, f.end
        note = "; ".join(f.notes) if f.notes else None
        src = "FOAM" if f.source.startswith("recovered") else "MITOS2"
        if f.type == "rep_origin":
            n += 1
            lines.append("\t".join([seqid, src, "origin_of_replication", str(s), str(e), ".",
                                    strand, ".", _gff_attr([("ID", "%s_ori%02d" % (prefix, n)),
                                                            ("Name", f.name), ("Note", note)])]))
            continue
        if f.type not in ("gene", "tRNA", "rRNA"):
            continue
        n += 1
        gid = "%s_g%02d" % (prefix, n)
        name = f.name + ("_%s" % chr(97 + f.part) if f.part is not None else "")
        if f.type == "gene":
            biotype, child = "protein_coding", "mRNA"
            product = PRODUCTS.get(f.name, f.name)
        elif f.type == "tRNA":
            biotype, child = "tRNA", "tRNA"
            aa = TRNA_AA.get(re.sub(r"^trn([A-Z]).*$", r"\1", f.name), "Xaa")
            product = "tRNA-" + aa
        else:
            biotype, child = "rRNA", "rRNA"
            product = PRODUCTS.get(f.name, f.name)
        base = [seqid, src, None, str(s), str(e), ".", strand]
        lines.append("\t".join(base[:2] + ["gene"] + base[3:] + ["."]) + "\t" +
                     _gff_attr([("ID", gid), ("Name", name), ("gene", f.name),
                                ("gene_biotype", biotype), ("Note", note)]))
        lines.append("\t".join(base[:2] + [child] + base[3:] + ["."]) + "\t" +
                     _gff_attr([("ID", gid + ".1"), ("Parent", gid), ("gene", f.name),
                                ("product", product),
                                ("anticodon", f.anticodon.lower() if f.anticodon else None)]))
        lines.append("\t".join(base[:2] + ["exon"] + base[3:] + ["."]) + "\t" +
                     _gff_attr([("ID", gid + ".1.exon1"), ("Parent", gid + ".1")]))
        if f.type == "gene":
            lines.append("\t".join(base[:2] + ["CDS"] + base[3:] + ["0"]) + "\t" +
                         _gff_attr([("ID", gid + ".1.cds"), ("Parent", gid + ".1"),
                                    ("product", product), ("transl_table", transl_table)]))
    with open(path, "w") as out:
        out.write("\n".join(lines) + "\n")


def write_bed(features, seqid, path, length=None):
    """BED6. With ``length``, a feature crossing the origin of the circular
    sequence becomes two lines with the same name (end .. length, 0 .. rest),
    since BED coordinates must lie within the sequence."""
    rows = []
    for f in features:
        sc = 0 if f.score is None else min(1000, int(round(-math.log10(f.score) * 10))
                                            if f.method == "mitfi" and f.score > 0
                                            else int(min(1000, f.score)))
        parts = [(f.start, f.end)]
        if length and f.end > length:
            parts = [(f.start, length), (0, f.end - length)]
        for s, e in parts:
            rows.append((s, e, f.label(), max(sc, 0), "+" if f.strand > 0 else "-"))
    with open(path, "w") as out:
        for r in sorted(rows):
            out.write("%s\t%d\t%d\t%s\t%d\t%s\n" % ((seqid,) + r))


def write_table(features, path):
    with open(path, "w") as out:
        out.write("#name\ttype\tstart\tend\tstrand\tlength\tscore\tanticodon\tpart\tcopy\tsource\tnotes\n")
        for f in sorted(features, key=lambda x: x.start):
            out.write("\t".join(str(x) for x in [
                f.name, f.type, f.start + 1, f.end, "+" if f.strand > 0 else "-", f.length,
                "." if f.score is None else f.score, f.anticodon or ".",
                "." if f.part is None else f.part, "." if f.copy is None else f.copy,
                f.source, "; ".join(f.notes) or "."]) + "\n")


def dump_json(obj, path):
    with open(path, "w") as out:
        json.dump(obj, out, indent=2, sort_keys=True)
        out.write("\n")
