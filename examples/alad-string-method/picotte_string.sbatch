#!/bin/bash
#SBATCH --job-name=cfacv-alad-string
#SBATCH --partition=def-sm
#SBATCH --account=abramsfreeprj
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Functional test of the cfacv string method under NAMD 3 on picotte.
#
# Alanine dipeptide in water, ~840 atoms, 24 replica systems (12 images x 2,
# SMPARAMS(dual) = 1).  One rank per replica fits on a single 48-core node,
# so the shipped alad_stringmethod.conf needs no editing.
#
# def-sm / abramsfreeprj is the free tier: no charge, lower priority, 12 h cap.
# For a billed run use --partition=def --account=<abramsprj|abramshkprj|
# abramsbbprj>; which billed account is Cameron's decision, so ask rather
# than infer it.
#
# Submit from the example directory:
#     cd <cfacv>/examples/alad-string-method && sbatch picotte_string.sbatch
#
set -euo pipefail

module purge
module load abramsGrp-namd/3.0.3/mpi-gcc

SRC=$SLURM_SUBMIT_DIR

# Where the cfacv tree lives.  Submitting from inside the repository makes
# $SRC/../.. correct, but the example is usually staged into a work area, so
# pass it in explicitly:
#     sbatch --export=ALL,CFACV_BASEDIR=/ifs/groups/abramsGrp/$USER/cfacv \
#            picotte_string.sbatch
if [ -z "${CFACV_BASEDIR:-}" ]; then
  CFACV_BASEDIR=$(cd "$SRC/../.." && pwd)
fi
if [ ! -d "$CFACV_BASEDIR/tcl" ]; then
  echo "ERROR: CFACV_BASEDIR='$CFACV_BASEDIR' is not a cfacv tree." >&2
  exit 1
fi

WORK=/scratch/$USER/$SLURM_JOB_ID

# MPI launches must work on a persistent filesystem: $TMPDIR and
# $BEEGFS_TMPDIR are destroyed the moment a job step starts.
mkdir -p "$WORK"
cp -a "$SRC"/. "$WORK"/
cd "$WORK"

# copy results back even if we are killed on walltime
cleanup() {
  echo "staging results back to $SRC/results.$SLURM_JOB_ID"
  mkdir -p "$SRC/results.$SLURM_JOB_ID"
  cp -a "$WORK"/output "$WORK"/*.log "$WORK"/*.history \
        "$SRC/results.$SLURM_JOB_ID"/ 2>/dev/null || true
}
trap cleanup EXIT

echo "=== node $(hostname), $SLURM_NTASKS ranks, work $WORK"
echo "=== namd3: $(command -v namd3)"

# --- stage 1: build the initial string ----------------------------------
# Restrained-MD warmups, one system at a time, no replicas, so this stage
# uses the single-node build.
if [ ! -f restr.job0.0.inp ]; then
  module purge
  module load abramsGrp-namd/3.0.3/multicore-gcc
  ./mkinitimages.sh --cfacv-base-dir="$CFACV_BASEDIR" --num-procs=8
  module purge
  module load abramsGrp-namd/3.0.3/mpi-gcc
fi

# --- stage 2: evolve the string -----------------------------------------
./runstring.sh --cfacv-base-dir="$CFACV_BASEDIR" --num-procs="$SLURM_NTASKS"

echo "=== done"
