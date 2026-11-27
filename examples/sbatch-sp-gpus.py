#!/usr/bin/env python3
"""Snakemake --cluster submit script for the CNAG SLURM cluster, sending
every job to the sp_gpus partition.

Usage:  --cluster "python3 sbatch-sp-gpus.py {dependencies}"

A stop-gap for when the CPU partitions (all, research, research_long) are
down or draining: there is no walltime-based routing, any "queue" in the
spec file is ignored, and no --time is passed, since sp_gpus has no time
limit.  No GPUs are requested; jobs only use the partition's CPUs and
memory.

Memory (MB): the spec file's "mem", else the rule's resources.mem_mb, else
the FOAM per-rule default below.  (sbatch-cnag.py falls back to 100 MB,
which is too little for flye or DIAMOND.)
"""

import os
import subprocess
import sys

from snakemake.utils import read_job_properties

PARTITION = "sp_gpus"
ACCOUNT = "sa_assembly"

# MB, used when neither the spec file nor the rule sets memory.
DEFAULT_MEM = 8000
RULE_MEM = {
    "flye_assembly": 64000,
    "polish": 32000,
    "prefilter_reads": 16000,
    "select_reads": 16000,
    "annotate_draft": 16000,
    "annotate_final": 16000,
}


def from_cluster(job_properties, key, default=None):
    """Spec-file value for key, else the job's own value, else default."""
    cluster = job_properties.get("cluster", {})
    if key in cluster:
        return cluster[key]
    return job_properties.get(key, default)


jobscript = sys.argv[-1]
dependencies = [d for d in sys.argv[1:-1] if d.strip()]

job_properties = read_job_properties(jobscript)
job_properties.setdefault("cluster", {})

rule = job_properties.get("rule", "snakejob")

cpus = from_cluster(job_properties, "threads", "1")
name = job_properties["cluster"].get(
    "name", job_properties.get("params", {}).get("name", "foam_" + rule)
)
mem = job_properties["cluster"].get("mem")
if not mem:
    mem = job_properties.get("resources", {}).get("mem_mb")
if not mem or str(mem) == "<TBD>":
    mem = RULE_MEM.get(rule, DEFAULT_MEM)
constraint = from_cluster(job_properties, "constraint", "")
array = from_cluster(job_properties, "array", "")

# A rule's log: directive is a list, and a rule without one gives an empty
# list, not a missing key.
log = job_properties["cluster"].get("log", job_properties.get("log") or [])
if isinstance(log, str):
    log = [log]

if len(log) >= 2:
    out, error = log[0], log[1]
elif len(log) == 1:
    out = error = log[0]
elif array:
    out, error = "logs/slurm-%A_%a.out", "logs/slurm-%A_%a.err"
else:
    out = "logs/%s-%%j.out" % rule
    error = "logs/%s-%%j.err" % rule

# sbatch accepts the job but SLURM fails it immediately if the log directory
# does not exist yet.
for path in (out, error):
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)

cmdline = ["sbatch", "--parsable"]
if dependencies:
    cmdline.append("--dependency=afterok:" + ",".join(dependencies))
cmdline += [
    "--job-name=%s" % name,
    "--output=%s" % out,
    "--error=%s" % error,
    "--partition=%s" % PARTITION,
    "--cpus-per-task=%s" % cpus,
    "--mem=%s" % mem,
    "--account=%s" % ACCOUNT,
]
# Only pass these when set: "--constraint=" / "--array=" are rejected by sbatch.
if constraint:
    cmdline.append("--constraint=%s" % constraint)
if array:
    cmdline.append("--array=%s" % array)
cmdline.append(jobscript)

try:
    result = subprocess.run(
        cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
except OSError as err:
    sys.stderr.write("sbatch failed for rule %s: %s\n" % (rule, err))
    sys.exit(1)

if result.returncode != 0:
    sys.stderr.write(
        "sbatch rejected rule %s (exit %d)\n  %s\n%s"
        % (rule, result.returncode, " ".join(cmdline), result.stderr)
    )
    sys.exit(result.returncode)

# --parsable prints "jobid" or "jobid;cluster"; snakemake wants just the id.
jobid = result.stdout.strip().split(";")[0]
if not jobid:
    sys.stderr.write("sbatch returned no job id for rule %s\n%s" % (rule, result.stderr))
    sys.exit(1)

print(jobid)
