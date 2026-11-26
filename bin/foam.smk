"""FOAM: assemble/finish complete circular mitogenomes from long reads.

Starting points (config key ``mode``):

  reads    pass 1: seed (bait mitogenome, or conserved mito proteins when no
                   bait) -> candidate reads -> select reads -> flye -> draft
           pass 2: the sample's own draft is the bait -> candidate reads ->
                   select reads -> flye -> draft        (fishing.second_pass)
           -> flye polish with the last pass's reads -> final annotate
  flye     flye contig(s) + reads: the contigs are the (self) bait ->
           candidate reads -> select reads -> draft -> polish -> final annotate
  hifiasm  hifiasm contig(s) or a whole assembly -> candidates -> final
           annotate (no reads, no polishing)

Read selection (fishing.method):
  kmer       permissive prefilter, then keep reads made of k-mers at
             mitochondrial coverage (rejects NUMT-carrying reads by coverage)
  alignment  keep reads whose alignments cover most of the read (FOAM 0.8)

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
KMER = FISH.get("kmer", {})
FLYE = config.get("flye", {})
MITOS = config.get("mitos", {})
CAND = config.get("candidates", {})
TAND = config.get("tandem", {})
THREADS = config.get("threads", {})
KEEP = config.get("keep_intermediate", False)
METHOD = FISH.get("method", "kmer")
if METHOD not in ("kmer", "alignment"):
    raise ValueError("fishing.method must be kmer or alignment")

if MODE in ("reads", "flye") and not READS:
    raise ValueError("mode %s needs inputs.reads" % MODE)
if MODE in ("flye", "hifiasm") and not INP.get("contigs"):
    raise ValueError("mode %s needs inputs.contigs" % MODE)


def bait_fasta():
    if INP.get("bait"):
        return INP["bait"]
    if INP.get("taxid") and INP.get("ncbi_email"):
        return "results/bait/bait.fa"
    return None


SEED = FISH.get("seed", "auto")
if SEED == "auto":
    SEED = "bait" if bait_fasta() else "genes"
if MODE == "reads" and SEED == "bait" and not bait_fasta():
    raise ValueError("fishing.seed is 'bait' but there is no inputs.bait or inputs.taxid+ncbi_email")

# last fishing pass, whose draft gets polished
LAST = 2 if (MODE == "reads" and FISH.get("second_pass", True)) else 1
FINAL = "results/final/%s_MT" % TOLID

wildcard_constraints:
    p="[12]",


def maybe_temp(path):
    return path if KEEP else temp(path)


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


# ----------------------------------------------------------------- seeds
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


rule seed_proteins:
    """Conserved mitochondrial proteins from the MITOS reference set: a seed
    that needs no related mitogenome and ignores gene order/control region."""
    output:
        fa="results/seed/proteins.fa",
        db="results/seed/proteins.dmnd",
    params:
        refdir=MITOS.get("refdir", os.path.join(BASE, "utils", "refseq_dir")),
        refseq=MITOS.get("refseq", "refseq81m"),
        genes=FISH.get("seed_genes", "cox1,cox2,cox3,cob,nad1,nad4,nad5"),
    log:
        "logs/seed_proteins.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "{FOAM} seed-proteins --refdir {params.refdir} --refseq {params.refseq} "
        "--genes {params.genes} -o {output.fa} > {log} 2>&1; "
        "diamond makedb --in {output.fa} -d {output.db} >> {log} 2>&1"


def bait_source(wc):
    if wc.p == "2":
        return "results/pass1/draft_MT.fa"          # the sample's own draft
    if MODE == "flye":
        return "results/pass1/candidates.fa"        # the user's flye contigs
    return bait_fasta()


rule double_bait:
    input:
        bait_source,
    output:
        maybe_temp("results/pass{p}/bait.doubled.fa"),
    log:
        "logs/double_bait.pass{p}.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "{FOAM} double-bait -i {input} -o {output} > {log} 2>&1"


# ----------------------------------------------------------------- fishing
def prefilter_inputs(wc):
    d = {"reads": READS}
    if wc.p == "1" and MODE == "reads" and SEED == "genes":
        d["db"] = rules.seed_proteins.output.db
    else:
        d["bait"] = "results/pass%s/bait.doubled.fa" % wc.p
    return d


def prefilter_cmd(wc, input, output, threads):
    """Shell for the first, cheap read screen.  With fishing.method=kmer it
    is deliberately permissive; the k-mer step does the real selection."""
    if "db" in input.keys():
        loop = " ".join(
            "diamond blastx -d {db} -q {r} --query-gencode {c} -F 15 --range-culling -k 1 "
            "--evalue 1e-5 --outfmt 6 qseqid -p {t} --quiet >> {raw};".format(
                db=input.db, r=r, c=config["genetic_code"], t=threads, raw=output.ids + ".raw")
            for r in input.reads)
        return "rm -f {raw}; {loop} sort -u {raw} > {ids}; rm -f {raw}".format(
            raw=output.ids + ".raw", loop=loop, ids=output.ids)
    if METHOD == "kmer":
        filt = "--min-read-frac 0 --min-identity 0 --min-aln {} --min-len {}".format(
            FISH.get("prefilter_min_aln", 300), FISH.get("min_len", 1000))
    else:
        filt = "--min-read-frac {} --min-aln {} --min-identity {} --min-len {} --max-len {}".format(
            FISH.get("min_read_frac", 0.8), FISH.get("min_aln", 1000),
            FISH.get("min_identity", 0.7), FISH.get("min_len", 1000), FISH.get("max_len", 0))
    return ("minimap2 -x {preset} --secondary=no -t {t} {bait} {reads} "
            "| {foam} fish-filter {filt} > {ids}").format(
        preset=FISH.get("minimap2_preset", "map-ont"), t=threads, bait=input.bait,
        reads=" ".join(input.reads), foam=FOAM, filt=filt, ids=output.ids)


rule prefilter_reads:
    """Stream all reads once; write only IDs and the (small) candidate set."""
    input:
        unpack(prefilter_inputs),
    output:
        ids=temp("results/pass{p}/candidate_ids.txt"),
        reads=temp("results/pass{p}/candidate_reads.fastq.gz"),
    params:
        cmd=prefilter_cmd,
    threads: THREADS.get("fish", 16)
    log:
        "logs/prefilter_reads.pass{p}.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "({params.cmd}) 2> {log}; "
        "seqkit grep -j {threads} -f {output.ids} {input.reads} -o {output.reads} 2>> {log}; "
        "echo \"$(wc -l < {output.ids}) candidate reads\" >> {log}"


rule select_reads:
    """k-mer abundance selection (method kmer) or pass-through (alignment)."""
    input:
        "results/pass{p}/candidate_reads.fastq.gz",
    output:
        reads=maybe_temp("results/pass{p}/mito_reads.fastq.gz"),
        report="results/pass{p}/read_selection.json",
    params:
        ids="results/pass{p}/mito_read_ids.txt",
        kopts="-k {} --sampling {} --coverage-fraction {} --min-score-fraction {}".format(
            KMER.get("k", 21), KMER.get("sampling", 8), KMER.get("coverage_fraction", 0.2),
            KMER.get("min_score_fraction", 0.5)),
    threads: 2
    log:
        "logs/select_reads.pass{p}.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        ("{FOAM} kmer-select {input} -o {params.ids} --report {output.report} {params.kopts} "
         "> {log} 2>&1; seqkit grep -f {params.ids} {input} -o {output.reads} 2>> {log}; "
         "rm -f {params.ids}") if METHOD == "kmer" else
        ("cp {input} {output.reads}; echo '{{\"method\": \"alignment\"}}' > {output.report}")


# ----------------------------------------------------------------- assemble
rule flye_assembly:
    input:
        reads="results/pass{p}/mito_reads.fastq.gz",
    output:
        fasta="results/pass{p}/flye/assembly.fasta",
        info="results/pass{p}/flye/assembly_info.txt",
        gfa="results/pass{p}/flye/assembly_graph.gfa",
    params:
        flag=FLYE.get("read_flag", "--nano-hq"),
        opts=FLYE.get("options", "--meta"),
        tmp="results/pass{p}/flye/run",
    threads: THREADS.get("flye", 16)
    log:
        "logs/flye_assembly.pass{p}.log",
    conda:
        "../envs/flye.yaml"
    shell:
        "rm -rf {params.tmp}; "
        "flye {params.flag} {input.reads} -o {params.tmp} -t {threads} {params.opts} > {log} 2>&1; "
        "for f in assembly.fasta assembly_info.txt assembly_graph.gfa flye.log; do "
        "mv {params.tmp}/$f results/pass{wildcards.p}/flye/; done; "
        "rm -rf {params.tmp}"


def candidate_inputs(wc):
    d = {}
    if MODE == "reads":
        d["fasta"] = "results/pass%s/flye/assembly.fasta" % wc.p
        d["info"] = "results/pass%s/flye/assembly_info.txt" % wc.p
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
    if INP.get("contig_names") and MODE != "reads":
        a += ["--names"] + list(INP["contig_names"])
    return " ".join(a)


rule candidates:
    input:
        unpack(candidate_inputs),
    output:
        fasta="results/pass{p}/candidates.fa",
        json="results/pass{p}/candidates.fa.json",
    params:
        args=candidate_args,
    threads: 4
    log:
        "logs/candidates.pass{p}.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "{FOAM} candidates -i {input.fasta} -o {output.fasta} -t {threads} {params.args} "
        "> {log} 2>&1"


# ----------------------------------------------------------------- annotate / polish
rule annotate_draft:
    """Split tandem copies, annotate, keep the most complete copy, orient."""
    input:
        "results/pass{p}/candidates.fa",
    output:
        fa="results/pass{p}/draft_MT.fa",
        report="results/pass{p}/draft_MT.report.json",
    params:
        args=annotate_args(),
        prefix="results/pass{p}/draft_MT",
        seqid="%s_MT_draft" % TOLID,
        work="results/pass{p}/annotate_work",
    threads: THREADS.get("annotate", 4)
    log:
        "logs/annotate_draft.pass{p}.log",
    conda:
        "../envs/mitos.yaml"
    shell:
        "{FOAM} annotate -i {input} -o {params.prefix} --seqid {params.seqid} "
        "--workdir {params.work} -t {threads} {params.args} > {log} 2>&1"


rule polish:
    """flye polishing of the chosen copy with the last pass's reads.  The
    copy is rotated by half its length first so the linear ends (poorly
    polished) are away from the start gene; the final step re-orients."""
    input:
        target="results/pass%d/draft_MT.fa" % LAST,
        reads="results/pass%d/mito_reads.fastq.gz" % LAST,
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
        "results/pass1/candidates.fa" if MODE == "hifiasm" else "results/polish/polished.fa",
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
