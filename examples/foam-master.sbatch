#!/bin/bash
# Snakemake master job for FOAM on SLURM (Snakemake >= 9).
#
# With Snakemake >= 9 the snakemake process stays alive for the whole run:
# it submits each job once its inputs exist and polls SLURM for job status.
# This script runs that master as a small, long job of its own.
#
# Submit it from the run directory (results/ and logs/ are written there),
# with the snakemake driver environment active (sbatch passes on PATH):
#
#   conda activate snakemake-9-slurm
#   export FOAM=/path/to/FOAM
#   export FOAM_CONDA_PREFIX=/path/to/foam_envs    # optional, shared conda envs
#   sbatch $FOAM/examples/foam-master.sbatch xgXerMont.foam.yaml
#
# Arguments: CONFIG [PROFILE] [extra snakemake options...]
#   PROFILE is a directory or a name under $FOAM/examples/profiles
#   (default sm9-slurm-cnag).
#
# The master must outlive the longest step, so it needs a partition without
# a (short) time limit. Override the partition on the command line, e.g.
# when only sp_gpus is up:
#
#   sbatch -p sp_gpus $FOAM/examples/foam-master.sbatch xgXerMont.foam.yaml sm9-slurm-sp_gpus
#
# Re-submitting after the master died resumes the run: finished outputs are
# skipped and the profiles set rerun-incomplete.

#SBATCH --job-name=foam_master
#SBATCH --partition=research_long
#SBATCH --account=sa_assembly
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=foam_master.%j.log

set -euo pipefail

CONFIG=${1:?usage: sbatch foam-master.sbatch CONFIG [PROFILE] [snakemake options...]}
PROFILE=${2:-sm9-slurm-cnag}
shift $(( $# >= 2 ? 2 : 1 ))

: "${FOAM:?set FOAM to the FOAM checkout before submitting (export FOAM=/path/to/FOAM)}"
if [[ ! -d $PROFILE ]]; then
    PROFILE=$FOAM/examples/profiles/$PROFILE
fi
[[ -f $PROFILE/config.yaml ]] || { echo "no profile config: $PROFILE/config.yaml" >&2; exit 1; }
[[ -f $CONFIG ]] || { echo "no config file: $CONFIG" >&2; exit 1; }
command -v snakemake >/dev/null || {
    echo "snakemake not on PATH: activate the driver environment before sbatch" >&2; exit 1; }

ARGS=(-s "$FOAM/bin/foam.smk" --configfile "$CONFIG" --profile "$PROFILE")
if [[ -n ${FOAM_CONDA_PREFIX:-} ]]; then
    ARGS+=(--conda-prefix "$FOAM_CONDA_PREFIX")
fi

echo "[$(date)] master job ${SLURM_JOB_ID:-local} on $(hostname) in $PWD"
echo "snakemake ${ARGS[*]} $*"
snakemake "${ARGS[@]}" "$@"
echo "[$(date)] done"
