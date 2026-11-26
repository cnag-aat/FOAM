"""Command-line entry point: ``python bin/foam.py <command> ...``"""

import argparse
import logging
import os
import shutil
import sys
from multiprocessing.pool import ThreadPool

from . import __version__
from . import annotation as ann
from . import candidates as cand
from . import fishing, seqio, tandem
from . import mitos as mt

log = logging.getLogger("foam")


# ---------------------------------------------------------------- annotate
def apply_taxon(a):
    """Fill genetic code / anchor / expected genes from --taxon; explicit
    options win."""
    preset = ann.taxon_preset(a.taxon) if a.taxon else {}
    if a.genetic_code is None:
        if not preset:
            sys.exit("give --taxon or --genetic-code")
        a.genetic_code = preset["genetic_code"]
    if not a.anchor:
        a.anchor = preset.get("anchor", ann.DEFAULT_ANCHOR)
    if not a.expected_genes:
        a.expected_genes = preset.get("genes", "metazoa")


def _annotate_one(job):
    (cid, seq, wdir, a) = job
    fa = os.path.join(wdir, cid + ".fa")
    seqio.write_fasta([(cid, seq)], fa)
    mdir = os.path.join(wdir, cid + ".mitos")
    if a.reuse_mitos and os.path.exists(os.path.join(mdir, "result.mitos")):
        log.info("reusing existing MITOS run %s", mdir)
    else:
        try:
            mt.run_mitos(fa, mdir, a.genetic_code, a.refdir, a.refseq,
                         runmitos=a.runmitos, extra_opts=a.mitos_opts)
        except mt.MitosError as e:
            log.error(str(e))
            return cid, None
    exp = ann.expected_genes(a.expected_genes)
    feats, q, changes = mt.annotate_dir(
        mdir, seq, exp, a.genetic_code, a.refdir, a.refseq, a.cmsearch,
        a.max_overlap, a.max_evalue, keep_origins=a.keep_origins,
        rescue=not a.no_rescue, threads=a.threads)
    return cid, (feats, q, changes, mdir)


def cmd_annotate(a):
    """Split tandem copies, annotate unique copies, keep the most complete,
    orient on the anchor gene and write FASTA/GFF3/BED/TSV/JSON."""
    apply_taxon(a)
    os.makedirs(a.workdir, exist_ok=True)
    records = seqio.read_fasta(a.input)
    if not records:
        sys.exit("no sequences in %s" % a.input)
    linear = set(a.linear or [])
    jobs, origin, seen = [], {}, {}
    contig_info = []
    for name, _d, seq in records:
        seq = seq.upper()
        circular = name not in linear
        if a.no_tandem:
            copies, info = [(0, len(seq), seq)], {"n_copies": 1, "length": len(seq)}
        else:
            copies, info = tandem.split_copies(seq, circular, a.min_period,
                                               a.minimap2, a.threads)
        info["contig"] = name
        contig_info.append(info)
        for i, (s, e, cs) in enumerate(copies):
            cid = "%s_copy%d" % (name, i + 1)
            key = cs
            if key in seen:            # identical copy: annotate once
                origin[cid] = (name, s, e, seen[key])
                continue
            seen[key] = cid
            origin[cid] = (name, s, e, cid)
            jobs.append((cid, cs, a.workdir, a))
    log.info("%d contig(s), %d unique copy/copies to annotate", len(records), len(jobs))
    pool = ThreadPool(max(1, min(a.parallel, len(jobs))))
    results = dict(pool.map(_annotate_one, jobs))
    pool.close()
    seqs = dict((j[0], j[1]) for j in jobs)

    table = []
    best = None
    for cid, (contig, s, e, rep) in sorted(origin.items()):
        r = results.get(rep)
        row = {"copy": cid, "contig": contig, "start": s + 1, "end": e,
               "length": e - s, "same_as": rep if rep != cid else None}
        if r is None:
            row["status"] = "mitos_failed"
        else:
            row["qc"] = r[1]
            if rep == cid and (best is None or ann.completeness_key(r[1]) >
                               ann.completeness_key(results[best][1])):
                best = cid
        table.append(row)
    if best is None:
        sys.exit("MITOS produced no usable annotation for any candidate")

    feats, q, changes, mdir = results[best]
    seq = seqs[best]
    oseq, ofeats, orient = mt.orient(seq, feats, a.anchor)
    if orient is None:
        log.warning("none of the anchor genes %s found; sequence left as is", a.anchor)
        orient = {"anchor": None}
    out = a.out_prefix
    d = os.path.dirname(out)
    if d:
        os.makedirs(d, exist_ok=True)
    seqid = a.seqid or os.path.basename(out)
    seqio.write_fasta([(seqid, "mitochondrion, complete genome, circular", oseq)], out + ".fa")
    prefix = a.gene_prefix or seqid.replace("_", "")
    ann.write_gff3(ofeats, seqid, len(oseq), out + ".gff3", prefix, a.genetic_code)
    ann.write_bed(ofeats, seqid, out + ".bed")
    ann.write_table(ofeats, out + ".genes.tsv")
    report = {"foam_version": __version__, "input": os.path.abspath(a.input),
              "contigs": contig_info, "candidates": table, "selected": best,
              "orientation": orient, "qc": q, "changes": changes}
    ann.dump_json(report, out + ".report.json")
    if a.keep_mitos:
        dest = out + ".mitos"
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(mdir, dest)
    if not a.keep_tmp:
        shutil.rmtree(a.workdir, ignore_errors=True)
    _print_summary(best, q, changes, orient)


def _print_summary(best, q, changes, orient):
    sys.stderr.write("selected %s (%d bp): %d/%d expected genes\n" % (
        best, q["length"], q["n_present"], q["n_expected"]))
    for k in ("missing", "fragmented", "duplicated", "recovered"):
        if q[k]:
            sys.stderr.write("  %s: %s\n" % (k, ", ".join(q[k])))
    for c in changes:
        sys.stderr.write("  fix: %s\n" % c)
    sys.stderr.write("  orientation: %s\n" % orient)


# ---------------------------------------------------------------- others
def cmd_split(a):
    recs = seqio.read_fasta(a.input)
    out = []
    for name, _d, seq in recs:
        copies, info = tandem.split_copies(seq.upper(), name not in (a.linear or []),
                                           a.min_period, a.minimap2, a.threads)
        sys.stderr.write("%s\n" % info)
        for i, (s, e, cs) in enumerate(copies):
            out.append(("%s_copy%d" % (name, i + 1), "%d-%d" % (s + 1, e), cs))
    seqio.write_fasta(out, a.output)


def cmd_rotate(a):
    recs = seqio.read_fasta(a.input)
    out = []
    for name, desc, seq in recs:
        off = len(seq) // 2 if a.half else a.offset
        out.append((name, desc, seqio.rotate(seq, off)))
    seqio.write_fasta(out, a.output)


def cmd_double_bait(a):
    seqio.write_fasta(fishing.double_bait(seqio.read_fasta(a.input)), a.output)


def cmd_fish_filter(a):
    src = sys.stdin if a.paf == "-" else seqio.open_text(a.paf)
    n = 0
    for rid in fishing.filter_paf(src, min_frac=a.min_read_frac, min_aln=a.min_aln,
                                  min_identity=a.min_identity, min_len=a.min_len,
                                  max_len=a.max_len):
        sys.stdout.write(rid + "\n")
        n += 1
    sys.stderr.write("fish-filter: %d reads kept\n" % n)


def cmd_candidates(a):
    recs = seqio.read_fasta(a.input)
    exp = a.expected_length
    if not exp and a.length_from and os.path.exists(a.length_from):
        lens = sorted(len(s) for _n, _d, s in seqio.read_fasta(a.length_from))
        exp = lens[len(lens) // 2] if lens else None
    if exp:
        if a.min_len is None:
            a.min_len = int(exp * (1 - a.length_tolerance))
        if a.max_len is None:
            a.max_len = int(exp * (1 + a.length_tolerance) * a.max_copies)
    a.min_len = a.min_len if a.min_len is not None else 5000
    a.max_len = a.max_len or 0
    circ = cand.flye_circular(a.flye_info) if a.flye_info else None
    cov = None
    if a.bait:
        cov = cand.bait_coverage(a.input, a.bait, a.minimap2, a.threads)
    keep, table = cand.select(recs, a.min_len, a.max_len, circ, not a.allow_linear,
                              a.names, cov, a.min_bait_cov)
    ann.dump_json(table, a.output + ".json")
    if not keep:
        sys.exit("no candidate contigs passed the filters; see %s.json" % a.output)
    seqio.write_fasta(keep, a.output)
    sys.stderr.write("%d candidate contig(s) written to %s\n" % (len(keep), a.output))


def cmd_fetch_bait(a):
    from . import ncbi
    res = ncbi.closest_mitogenomes(a.taxid, a.email, max_results=a.max_results,
                                   api_key=a.api_key)
    seqio.write_fasta([(r["accession"], r["title"], r["sequence"]) for r in res["records"]],
                      a.output)
    ann.dump_json(dict((k, v) for k, v in res.items() if k != "records"), a.output + ".json")


def build_parser():
    p = argparse.ArgumentParser(prog="foam", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("annotate", help="split tandem copies, annotate, pick best, orient")
    s.add_argument("-i", "--input", required=True, help="candidate contig(s) FASTA")
    s.add_argument("-o", "--out-prefix", required=True)
    s.add_argument("--taxon", help="taxon preset; sets genetic code, start gene and "
                   "expected genes: " + ", ".join(sorted(ann.TAXON_PRESETS)))
    s.add_argument("-c", "--genetic-code", type=int,
                   help="mitochondrial genetic code (overrides --taxon)")
    s.add_argument("--refdir", required=True, help="MITOS reference directory (parent)")
    s.add_argument("--refseq", default="refseq81m")
    s.add_argument("--runmitos", default="runmitos.py")
    s.add_argument("--mitos-opts", default="", help="extra runmitos options")
    s.add_argument("--anchor", nargs="+",
                   help="gene(s) whose 5' end becomes position 1, in order of preference "
                        "(overrides --taxon; default cox1)")
    s.add_argument("--expected-genes", help="gene set name or comma list (overrides --taxon)")
    s.add_argument("--seqid", help="output sequence name (default: basename of prefix)")
    s.add_argument("--gene-prefix")
    s.add_argument("--linear", nargs="*", help="contig names that are NOT circular")
    s.add_argument("--no-tandem", action="store_true")
    s.add_argument("--min-period", type=int, default=5000)
    s.add_argument("--max-overlap", type=int, default=15)
    s.add_argument("--max-evalue", type=float, default=0.01)
    s.add_argument("--minimap2", default="minimap2")
    s.add_argument("--workdir", default="foam_annotate_tmp")
    s.add_argument("--parallel", type=int, default=2, help="MITOS runs in parallel")
    s.add_argument("-t", "--threads", type=int, default=2)
    s.add_argument("--cmsearch", default="cmsearch", help="for rescuing missed tRNAs")
    s.add_argument("--no-rescue", action="store_true",
                   help="do not recover/rescue genes MITOS dropped")
    s.add_argument("--keep-origins", action="store_true",
                   help="keep MITOS OH/OL predictions in the output")
    s.add_argument("--reuse-mitos", action="store_true",
                   help="reuse MITOS results already present in --workdir")
    s.add_argument("--keep-mitos", action="store_true", help="copy the selected MITOS run")
    s.add_argument("--keep-tmp", action="store_true")
    s.set_defaults(func=cmd_annotate)

    s = sub.add_parser("split-copies", help="write tandem copies of each contig")
    s.add_argument("-i", "--input", required=True)
    s.add_argument("-o", "--output", required=True)
    s.add_argument("--linear", nargs="*")
    s.add_argument("--min-period", type=int, default=5000)
    s.add_argument("--minimap2", default="minimap2")
    s.add_argument("-t", "--threads", type=int, default=2)
    s.set_defaults(func=cmd_split)

    s = sub.add_parser("rotate", help="rotate circular sequences")
    s.add_argument("-i", "--input", required=True)
    s.add_argument("-o", "--output", required=True)
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--offset", type=int)
    g.add_argument("--half", action="store_true")
    s.set_defaults(func=cmd_rotate)

    s = sub.add_parser("double-bait", help="concatenate each bait with itself")
    s.add_argument("-i", "--input", required=True)
    s.add_argument("-o", "--output", required=True)
    s.set_defaults(func=cmd_double_bait)

    s = sub.add_parser("fish-filter", help="PAF (stdin) -> mitochondrial read IDs")
    s.add_argument("paf", nargs="?", default="-")
    s.add_argument("--min-read-frac", type=float, default=0.8)
    s.add_argument("--min-aln", type=int, default=1000)
    s.add_argument("--min-identity", type=float, default=0.7)
    s.add_argument("--min-len", type=int, default=1000)
    s.add_argument("--max-len", type=int, default=0)
    s.set_defaults(func=cmd_fish_filter)

    s = sub.add_parser("candidates", help="select candidate mito contigs")
    s.add_argument("-i", "--input", required=True)
    s.add_argument("-o", "--output", required=True)
    s.add_argument("--flye-info")
    s.add_argument("--names", nargs="*")
    s.add_argument("--expected-length", type=int)
    s.add_argument("--length-from", help="FASTA (e.g. bait) whose median length is expected")
    s.add_argument("--length-tolerance", type=float, default=0.25)
    s.add_argument("--max-copies", type=int, default=4,
                   help="allow contigs up to this many tandem copies long")
    s.add_argument("--min-len", type=int)
    s.add_argument("--max-len", type=int)
    s.add_argument("--allow-linear", action="store_true")
    s.add_argument("--bait")
    s.add_argument("--min-bait-cov", type=float, default=0.3)
    s.add_argument("--minimap2", default="minimap2")
    s.add_argument("-t", "--threads", type=int, default=4)
    s.set_defaults(func=cmd_candidates)

    s = sub.add_parser("fetch-bait", help="closest RefSeq/GenBank mitogenome(s) by taxid")
    s.add_argument("--taxid", required=True)
    s.add_argument("--email", required=True)
    s.add_argument("--api-key")
    s.add_argument("--max-results", type=int, default=3)
    s.add_argument("-o", "--output", required=True)
    s.set_defaults(func=cmd_fetch_bait)
    return p


def main(argv=None):
    p = build_parser()
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="[foam %(levelname)s] %(message)s")
    if not getattr(a, "func", None):
        p.print_help()
        return 1
    a.func(a)
    return 0
