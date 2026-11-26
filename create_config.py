#!/usr/bin/env python3
"""Write a FOAM config file.

Three starting points (``--mode``):

  reads     reads -> fish mito reads (bait mitogenome or conserved proteins)
            -> flye -> re-fish with own draft -> flye -> polish -> annotate
  flye      flye circular contig(s) + reads   -> polish -> annotate
  hifiasm   hifiasm contig(s) (or a whole assembly) -> annotate (no polishing,
            no reads needed)

Taxon presets set the genetic code, the start gene (trnF for vertebrates, the
5' end of cox1 for everything else) and the expected gene set.  With
``--taxid`` the preset and genetic code are looked up in NCBI Taxonomy.

Examples
--------
  create_config.py --mode hifiasm --tolid xgXerMont --taxon mollusca \\
      --contigs mito_candidate.fa -o xgXerMont.foam.yaml

  create_config.py --mode reads --tolid ibLobCana --taxid 1234 --email me@x.org \\
      --reads reads.ont.fastq.gz --read-type ont-hq -o ibLobCana.foam.yaml
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from foam import annotation as ann  # noqa: E402

READ_TYPES = {
    # flye input flag, minimap2 preset for fishing
    "ont-hq": ("--nano-hq", "map-ont"),
    "ont-raw": ("--nano-raw", "map-ont"),
    "ont-corr": ("--nano-corr", "map-ont"),
    "hifi": ("--pacbio-hifi", "map-hifi"),
    "pacbio-raw": ("--pacbio-raw", "map-pb"),
}


def yaml_dump(obj, indent=0):
    """Tiny YAML emitter for nested dicts/lists of scalars (no PyYAML needed)."""
    pad = "  " * indent
    lines = []
    for k, v in obj.items():
        if isinstance(v, dict):
            lines.append("%s%s:" % (pad, k))
            lines.append(yaml_dump(v, indent + 1))
        elif isinstance(v, list):
            lines.append("%s%s: %s" % (pad, k, json.dumps(v)))
        else:
            lines.append("%s%s: %s" % (pad, k, scalar(v)))
    return "\n".join(lines)


def scalar(v):
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    return json.dumps(str(v))


def abspath_or_die(p, what):
    if p is None:
        return None
    ap = os.path.abspath(os.path.expanduser(p))
    if not os.path.exists(ap):
        sys.exit("error: %s not found: %s" % (what, p))
    return ap


def build_config(a):
    preset_name = a.taxon
    code = a.genetic_code
    if a.taxid and a.email and not a.offline and (preset_name is None or code is None):
        from foam import ncbi
        tax = ncbi.taxonomy(a.taxid, ncbi.Entrez(a.email, a.api_key))
        preset_name = preset_name or ncbi.preset_for_lineage(tax["lineage"])
        code = code or tax["mito_genetic_code"]
        sys.stderr.write("NCBI: %s -> preset %s, mito genetic code %s\n" % (
            tax["taxon"]["ScientificName"], preset_name, code))
    if preset_name is None:
        sys.exit("error: give --taxon (or --taxid with --email)")
    preset = ann.taxon_preset(preset_name)

    reads = [abspath_or_die(r, "reads") for r in (a.reads or [])]
    contigs = abspath_or_die(a.contigs, "contigs")
    if a.mode in ("reads", "flye") and not reads:
        sys.exit("error: --mode %s needs --reads" % a.mode)
    if a.mode in ("flye", "hifiasm") and not contigs:
        sys.exit("error: --mode %s needs --contigs" % a.mode)
    if a.mode == "reads" and a.seed == "bait" and not (a.bait or (a.taxid and a.email)):
        sys.exit("error: --seed bait needs --bait, or --taxid with --email to fetch one")
    if a.taxid and not a.bait and a.mode == "reads" and not a.email:
        sys.exit("error: fetching a bait from NCBI needs --email")
    flye_flag, mm2_preset = READ_TYPES[a.read_type]

    cfg = {
        "tolid": a.tolid,
        "mode": a.mode,
        "taxon": preset_name,
        "genetic_code": code if code is not None else preset["genetic_code"],
        "anchor": a.anchor or preset["anchor"],
        "expected_genes": a.expected_genes or preset["genes"],
        "expected_length": a.expected_length,
        "inputs": {
            "reads": reads,
            "read_type": a.read_type,
            "contigs": contigs,
            "contig_names": a.contig_names or [],
            "flye_info": abspath_or_die(a.flye_info, "flye assembly_info"),
            "bait": abspath_or_die(a.bait, "bait"),
            "taxid": a.taxid,
            "ncbi_email": a.email,
        },
        "candidates": {
            "require_circular": not a.allow_linear,
            "max_copies": 4,
            "length_tolerance": 0.25,
            "min_bait_cov": 0.3,
        },
        "fishing": {
            # seed for the first pass: bait (a related mitogenome), genes
            # (conserved mito proteins, no bait needed) or auto
            "seed": a.seed,
            "seed_genes": "cox1,cox2,cox3,cob,nad1,nad4,nad5",
            # kmer: permissive prefilter, then reads made of mito-coverage
            # k-mers; alignment: reads mostly covered by bait alignments
            "method": a.fishing_method,
            "second_pass": not a.no_second_pass,
            "prefilter_min_aln": 300,
            # backend: fastk (Gene Myers' FastK; all k-mers, fast), python
            # (1/sampling of the k-mers, no dependencies) or auto
            "kmer": {"backend": "auto", "k": 21, "sampling": 8, "coverage_fraction": 0.2,
                     "min_score_fraction": 0.5, "memory_gb": 4},
            "minimap2_preset": mm2_preset,
            "min_read_frac": 0.8,
            "min_aln": 1000,
            "min_identity": 0.7,
            "min_len": 1000,
            "max_len": 0,
        },
        "flye": {"read_flag": flye_flag, "options": "--meta", "polish_iterations": 1},
        "tandem": {"resolve": True, "min_period": 5000},
        "mitos": {
            "refdir": os.path.abspath(a.refdir),
            "refseq": a.refseq,
            "runmitos": a.runmitos,
            "options": "",
            "max_overlap": 15,
            "max_evalue": 0.01,
            "parallel": 2,
        },
        "threads": {"fish": a.threads, "kmer": min(a.threads, 8), "flye": a.threads,
                    "annotate": 4},
        "keep_intermediate": False,
    }
    return cfg


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", required=True, choices=["reads", "flye", "hifiasm"])
    p.add_argument("--tolid", required=True)
    p.add_argument("-o", "--output", help="config file (default <tolid>.foam.yaml)")
    g = p.add_argument_group("taxon")
    g.add_argument("--taxon", choices=sorted(ann.TAXON_PRESETS))
    g.add_argument("--taxid", help="NCBI taxid (preset/genetic code lookup and bait fetch)")
    g.add_argument("--email", help="e-mail for NCBI E-utilities")
    g.add_argument("--api-key")
    g.add_argument("--offline", action="store_true", help="do not query NCBI now")
    g.add_argument("--genetic-code", type=int, help="override the preset genetic code")
    g.add_argument("--anchor", nargs="+", help="override the start gene(s)")
    g.add_argument("--expected-genes", help="gene set name or comma-separated list")
    g.add_argument("--expected-length", type=int, help="approx. mitogenome length")
    g = p.add_argument_group("inputs")
    g.add_argument("--reads", nargs="+")
    g.add_argument("--read-type", default="ont-hq", choices=sorted(READ_TYPES))
    g.add_argument("--contigs", help="candidate contig(s) or whole assembly FASTA")
    g.add_argument("--contig-names", nargs="+", help="use only these contigs")
    g.add_argument("--flye-info", help="flye assembly_info.txt (circularity)")
    g.add_argument("--allow-linear", action="store_true")
    g.add_argument("--bait", help="bait/seed mitogenome FASTA")
    g = p.add_argument_group("read fishing (reads/flye modes)")
    g.add_argument("--seed", default="auto", choices=["auto", "bait", "genes"],
                   help="first-pass seed; auto = bait if given/fetchable, else genes")
    g.add_argument("--fishing-method", default="kmer", choices=["kmer", "alignment"])
    g.add_argument("--no-second-pass", action="store_true",
                   help="do not re-fish with the sample's own draft (reads mode)")
    g = p.add_argument_group("software")
    g.add_argument("--refdir", default=os.path.join(HERE, "utils", "refseq_dir"))
    g.add_argument("--refseq", default="refseq81m")
    g.add_argument("--runmitos", default="runmitos.py",
                   help="MITOS command (runmitos.py for 2.1.3, runmitos for >=2.1.7)")
    g.add_argument("--threads", type=int, default=16)
    a = p.parse_args(argv)
    cfg = build_config(a)
    out = a.output or "%s.foam.yaml" % a.tolid
    with open(out, "w") as fh:
        fh.write("# FOAM config written by create_config.py\n")
        fh.write(yaml_dump(cfg) + "\n")
    sys.stderr.write("wrote %s\n" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
