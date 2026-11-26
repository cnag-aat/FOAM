"""FOAM: assemble/finish complete circular mitogenomes from long reads.

Starting points (config key ``mode``):

  reads    bait -> fish mito reads -> flye -> candidate contigs -> annotate
           (pick best tandem copy, orient) -> flye polish -> final annotate
  flye     flye contig(s) + reads   -> candidates -> fish reads with the
           contigs as bait -> annotate -> polish -> final annotate
  hifiasm  hifiasm contig(s) or a whole hifiasm assembly -> candidates ->
           final annotate (no polishing, no reads)

Write a config with ``create_config.py`` and run e.g.

  snakemake -s bin/foam.smk --configfile sp.foam.yaml --use-conda -j 16
"""

import os

BASE = os.path.abspath(os.path.join(workflow.basedir, ".."))
FOAM = "python {}/bin/foam.py".format(BASE)

MODE = config["mode"]
if MODE not in ("reads", "flye", "hifiasm"):
    raise ValueError("config 'mode' must be reads, flye or hifiasm (got %r)" % MODE)

TOLID = config["tolid"]
INP = config.get("inputs", {})
READS = INP.get("reads") or []
if isinstance(READS, str):
    READS = [READS]
FISH = config.get("fishing", {})
FLYE = config.get("flye", {})
MITOS = config.get("mitos", {})
CAND = config.get("candidates", {})
TAND = config.get("tandem", {})
THREADS = config.get("threads", {})
KEEP = config.get("keep_intermediate", False)

if MODE in ("reads", "flye") and not READS:
    raise ValueError("mode %s needs inputs.reads" % MODE)
if MODE in ("flye", "hifiasm") and not INP.get("contigs"):
    raise ValueError("mode %s needs inputs.contigs" % MODE)

FINAL = "results/final/%s_MT" % TOLID


def maybe_temp(path):
    return path if KEEP else temp(path)


def bait_fasta():
    if INP.get("bait"):
        return INP["bait"]
    if INP.get("taxid") and INP.get("ncbi_email"):
        return "results/bait/bait.fa"
    return None


def annotate_args():
    a = ["--taxon", config["taxon"],
         "-c", str(config["genetic_code"]),
         "--anchor"] + list(config["anchor"]) + [
         "--expected-genes", str(config.get("expected_genes", "metazoa")),
         "--refdir", MITOS.get("refdir", os.path.join(BASE, "utils", "refseq_dir")),
         "--refseq", MITOS.get("refseq", "refseq81m"),
         "--runmitos", "'%s'" % MITOS.get("runmitos", "runmitos.py"),
         "--max-overlap", str(MITOS.get("max_overlap", 15)),
         "--max-evalue", str(MITOS.get("max_evalue", 0.01)),
         "--parallel", str(MITOS.get("parallel", 2)),
         "--min-period", str(TAND.get("min_period", 5000))]
    if MITOS.get("options"):
        a += ["--mitos-opts", "'%s'" % MITOS["options"]]
    if not TAND.get("resolve", True):
        a.append("--no-tandem")
    if MITOS.get("keep_origins"):
        a.append("--keep-origins")
    return " ".join(a)


rule all:
    input:
        expand(FINAL + ".{ext}", ext=["fa", "gff3", "bed", "genes.tsv", "report.json"]),


localrules: all, fetch_bait


# ----------------------------------------------------------------- bait
rule fetch_bait:
    """Closest complete mitogenome(s) in NCBI (needs internet: local rule)."""
    output:
        "results/bait/bait.fa",
    params:
        taxid=INP.get("taxid"),
        email=INP.get("ncbi_email"),
    log:
        "logs/fetch_bait.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "{FOAM} fetch-bait --taxid {params.taxid} --email {params.email} "
        "-o {output} > {log} 2>&1"


def fishing_bait(wc):
    if MODE == "flye":
        return "results/candidates/candidates.fa"
    b = bait_fasta()
    if b is None:
        raise ValueError("mode reads needs inputs.bait or inputs.taxid + inputs.ncbi_email")
    return b


rule double_bait:
    input:
        fishing_bait,
    output:
        maybe_temp("results/fishing/bait.doubled.fa"),
    log:
        "logs/double_bait.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "{FOAM} double-bait -i {input} -o {output} > {log} 2>&1"


rule fish_reads:
    """Stream reads through minimap2 against the doubled bait; keep reads whose
    alignments cover most of the read.  Nothing but read IDs and the (small)
    selected read set is written."""
    input:
        bait=rules.double_bait.output[0],
        reads=READS,
    output:
        reads=maybe_temp("results/fishing/mito_reads.fastq.gz"),
        ids=temp("results/fishing/mito_read_ids.txt"),
    params:
        preset=FISH.get("minimap2_preset", "map-ont"),
        filt="--min-read-frac {} --min-aln {} --min-identity {} --min-len {} --max-len {}".format(
            FISH.get("min_read_frac", 0.8), FISH.get("min_aln", 1000),
            FISH.get("min_identity", 0.7), FISH.get("min_len", 1000), FISH.get("max_len", 0)),
    threads: THREADS.get("fish", 16)
    log:
        "logs/fish_reads.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "(minimap2 -x {params.preset} --secondary=no -t {threads} {input.bait} {input.reads} "
        "| {FOAM} fish-filter {params.filt} > {output.ids}) 2> {log}; "
        "seqkit grep -j {threads} -f {output.ids} {input.reads} -o {output.reads} 2>> {log}; "
        "echo \"$(wc -l < {output.ids}) reads kept\" >> {log}"


# ----------------------------------------------------------------- assemble
rule flye_assembly:
    input:
        reads="results/fishing/mito_reads.fastq.gz",
    output:
        fasta="results/flye/assembly.fasta",
        info="results/flye/assembly_info.txt",
        gfa="results/flye/assembly_graph.gfa",
    params:
        flag=FLYE.get("read_flag", "--nano-hq"),
        opts=FLYE.get("options", "--meta"),
        tmp="results/flye/run",
    threads: THREADS.get("flye", 16)
    log:
        "logs/flye_assembly.log",
    conda:
        "../envs/flye.yaml"
    shell:
        "rm -rf {params.tmp}; "
        "flye {params.flag} {input.reads} -o {params.tmp} -t {threads} {params.opts} > {log} 2>&1; "
        "mv {params.tmp}/assembly.fasta {output.fasta}; "
        "mv {params.tmp}/assembly_info.txt {output.info}; "
        "mv {params.tmp}/assembly_graph.gfa {output.gfa}; "
        "mv {params.tmp}/flye.log results/flye/flye.log; "
        "rm -rf {params.tmp}"


def candidate_inputs(wc):
    d = {}
    if MODE == "reads":
        d["fasta"] = rules.flye_assembly.output.fasta
        d["info"] = rules.flye_assembly.output.info
    else:
        d["fasta"] = INP["contigs"]
        if INP.get("flye_info"):
            d["info"] = INP["flye_info"]
    b = bait_fasta()
    if b and MODE != "flye":
        d["bait"] = b
    return d


def candidate_args(wc, input):
    a = ["--length-tolerance", str(CAND.get("length_tolerance", 0.25)),
         "--max-copies", str(CAND.get("max_copies", 4)),
         "--min-bait-cov", str(CAND.get("min_bait_cov", 0.3))]
    if config.get("expected_length"):
        a += ["--expected-length", str(config["expected_length"])]
    if "bait" in input.keys():
        a += ["--length-from", input.bait]
        if MODE == "hifiasm":           # screen a whole assembly against the bait
            a += ["--bait", input.bait]
    if "info" in input.keys():
        a += ["--flye-info", input.info]
    elif MODE == "flye" or not CAND.get("require_circular", True):
        a.append("--allow-linear")
    if INP.get("contig_names"):
        a += ["--names"] + list(INP["contig_names"])
    return " ".join(a)


rule candidates:
    input:
        unpack(candidate_inputs),
    output:
        fasta="results/candidates/candidates.fa",
        json="results/candidates/candidates.fa.json",
    params:
        args=candidate_args,
    threads: 4
    log:
        "logs/candidates.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "{FOAM} candidates -i {input.fasta} -o {output.fasta} -t {threads} {params.args} "
        "> {log} 2>&1"


# ----------------------------------------------------------------- annotate / polish
rule annotate_prepolish:
    """Split tandem copies, annotate, keep the most complete copy, orient."""
    input:
        "results/candidates/candidates.fa",
    output:
        fa="results/prepolish/%s_MT.fa" % TOLID,
        report="results/prepolish/%s_MT.report.json" % TOLID,
    params:
        args=annotate_args(),
        prefix="results/prepolish/%s_MT" % TOLID,
        work="results/prepolish/work",
    threads: THREADS.get("annotate", 4)
    log:
        "logs/annotate_prepolish.log",
    conda:
        "../envs/mitos.yaml"
    shell:
        "{FOAM} annotate -i {input} -o {params.prefix} --workdir {params.work} "
        "-t {threads} {params.args} > {log} 2>&1"


rule polish:
    """flye polishing of the single chosen copy.  The copy is rotated by half
    its length first so the linear ends (poorly polished) are not at the
    start gene; the final step re-orients."""
    input:
        target=rules.annotate_prepolish.output.fa,
        reads="results/fishing/mito_reads.fastq.gz",
    output:
        "results/polish/polished.fa",
    params:
        flag=FLYE.get("read_flag", "--nano-hq"),
        iters=FLYE.get("polish_iterations", 1),
        tmp="results/polish/run",
        rotated="results/polish/target.rotated.fa",
    threads: THREADS.get("flye", 16)
    log:
        "logs/polish.log",
    conda:
        "../envs/flye.yaml"
    shell:
        "rm -rf {params.tmp}; "
        "python {BASE}/bin/foam.py rotate --half -i {input.target} -o {params.rotated} 2> {log}; "
        "flye --polish-target {params.rotated} {params.flag} {input.reads} "
        "-i {params.iters} -o {params.tmp} -t {threads} >> {log} 2>&1; "
        "mv {params.tmp}/polished_{params.iters}.fasta {output}; "
        "rm -rf {params.tmp} {params.rotated}"


rule annotate_final:
    input:
        "results/candidates/candidates.fa" if MODE == "hifiasm" else "results/polish/polished.fa",
    output:
        expand(FINAL + ".{ext}", ext=["fa", "gff3", "bed", "genes.tsv", "report.json"]),
    params:
        args=annotate_args(),
        prefix=FINAL,
        seqid="%s_MT" % TOLID,
        work="results/final/work",
        keep="--keep-mitos" if KEEP else "",
    threads: THREADS.get("annotate", 4)
    log:
        "logs/annotate_final.log",
    conda:
        "../envs/mitos.yaml"
    shell:
        "{FOAM} annotate -i {input} -o {params.prefix} --seqid {params.seqid} "
        "--workdir {params.work} -t {threads} {params.keep} {params.args} > {log} 2>&1"
