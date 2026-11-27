#!/usr/bin/env python3
"""Snakemake --cluster submit script for the CNAG SLURM cluster.

Usage:  --cluster "python3 sbatch-cnag.py {dependencies}"

Partition is chosen from the requested walltime, since each partition
enforces its own limit and sbatch rejects jobs that ask for more:

    <= 1:00:00    -> all            (partition limit 1:00:00)
    <= 12:00:00   -> research       (partition limit 12:00:00)
    >  12:00:00   -> research_long  (no limit)

A "queue" in the spec file naming some other partition (e.g. a GPU one)
is honoured as-is; the rule above only knows about the CPU partitions.
"""

import os
import re
import subprocess
import sys

from snakemake.utils import read_job_properties

# Partitions selected by walltime, longest-fitting first.
TIME_ROUTED_PARTITIONS = ("all", "research", "research_long")
PARTITION_LIMITS = (
    (3600, "all"),           #  1 h
    (43200, "research"),     # 12 h
)
LONG_PARTITION = "research_long"


def parse_time(spec):
    """Return walltime in seconds for any SLURM time format.

    Accepts minutes, minutes:seconds, hours:minutes:seconds,
    days-hours, days-hours:minutes, days-hours:minutes:seconds,
    and infinite/unlimited.
    """
    spec = str(spec).strip()
    if spec.lower() in ("infinite", "unlimited"):
        return None  # no limit -> always the long partition

    days = 0
    if "-" in spec:
        d, _, spec = spec.partition("-")
        days = int(d)

    parts = spec.split(":") if spec else ["0"]
    if not all(re.fullmatch(r"\d+", p or "0") for p in parts):
        raise ValueError("unparsable time: %r" % spec)
    parts = [int(p or 0) for p in parts]

    if days:
        # days-hours[:minutes[:seconds]]
        h, m, s = (parts + [0, 0, 0])[:3]
    elif len(parts) == 1:
        h, m, s = 0, parts[0], 0            # bare number == minutes
    elif len(parts) == 2:
        h, m, s = 0, parts[0], parts[1]     # minutes:seconds
    elif len(parts) == 3:
        h, m, s = parts                     # hours:minutes:seconds
    else:
        raise ValueError("unparsable time: %r" % spec)

    return days * 86400 + h * 3600 + m * 60 + s


def choose_partition(time_spec, requested):
    """Pick the partition for this walltime, honouring non-CPU queues."""
    if requested and requested not in TIME_ROUTED_PARTITIONS:
        return requested

    try:
        seconds = parse_time(time_spec)
    except (ValueError, TypeError):
        return LONG_PARTITION  # can't tell -> the partition with no limit

    if seconds is None:
        return LONG_PARTITION
    for limit, partition in PARTITION_LIMITS:
        if seconds <= limit:
            return partition
    return LONG_PARTITION


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
    "name", job_properties.get("params", {}).get("name", "Snakejob")
)
mem = from_cluster(job_properties, "mem", "100")
constraint = from_cluster(job_properties, "constraint", "")
array = from_cluster(job_properties, "array", "")

time = job_properties["cluster"].get(
    "time", job_properties.get("params", {}).get("time", "10:00:00")
)
partition = choose_partition(time, from_cluster(job_properties, "queue", ""))

# Log files. A rule's log: directive is a list, and most CLAWS rules have no
# log: at all -- which makes job_properties["log"] an empty list, not a
# missing key. Indexing it blindly is what stopped those rules submitting.
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
    "--time=%s" % time,
    "--partition=%s" % partition,
    "--cpus-per-task=%s" % cpus,
    "--mem=%s" % mem,
    "--account=sa_assembly",
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
