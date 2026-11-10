#!/bin/bash
#
# String method in collective variables via NAMD/replica
#
# (c) 2016 Cameron F Abrams, Drexel University
#
# Reworked 2026-11-10 to resolve NAMD and cfacv from the environment instead
# of hardcoded paths on the original author's workstation, and to launch the
# NAMD 3 MPI build, which is what provides +replicas across nodes.
#
# Needs an MPI (or other multi-partition) NAMD build: a multicore build
# refuses +replicas with "+partitions other than 1 is not allowed for
# multicore build".  With namd3 and mpirun on PATH:
#
#     ./runstring.sh --num-procs=<total ranks>
#
# NAMD's +replicas requires the total process count to be an exact multiple
# of the replica count; this script checks that before launching.
#
set -euo pipefail

CFACV_BASEDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
NAMD3=namd3
MPIRUN="mpirun"
MPIRUN_OPTS="--map-by core --bind-to core"
CONF="job0.conf"
STRING_METHOD_CONFIG="alad_stringmethod.conf"
NP=""          # total MPI ranks; default = num_replicas (one rank per replica)
NI=""          # number of replicas; default = num_replicas from the config
OUTDIR="output"

for i in "$@"; do
  case $i in
    -np=*|--num-procs=*)            NP="${i#*=}" ;;
    -ni=*|--num-replicas=*)         NI="${i#*=}" ;;
    -c=*|--conf=*)                  CONF="${i#*=}" ;;
    -i=*|--string-method-config=*)  STRING_METHOD_CONFIG="${i#*=}" ;;
    -cfacvbd=*|--cfacv-base-dir=*)  CFACV_BASEDIR="${i#*=}" ;;
    -namd3=*|--namd3-command=*)     NAMD3="${i#*=}" ;;
    -mpirun=*|--mpirun-command=*)   MPIRUN="${i#*=}" ;;
    --mpirun-opts=*)                MPIRUN_OPTS="${i#*=}" ;;
    -o=*|--output-dir=*)            OUTDIR="${i#*=}" ;;
    -h|--help)
      sed -n '2,20p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "runstring.sh: unknown option '$i'" >&2; exit 2 ;;
  esac
done

for f in "$CONF" "$STRING_METHOD_CONFIG"; do
  [ -f "$f" ] || { echo "ERROR: file $f not found." >&2; exit 1; }
done
[ -d "$CFACV_BASEDIR/tcl" ] || { echo "ERROR: $CFACV_BASEDIR is not a cfacv tree." >&2; exit 1; }
command -v "$NAMD3"  >/dev/null || { echo "ERROR: '$NAMD3' not on PATH. module load abramsGrp-namd/3.0.3/mpi-gcc" >&2; exit 1; }
command -v "$MPIRUN" >/dev/null || { echo "ERROR: '$MPIRUN' not on PATH." >&2; exit 1; }

# replica count comes from the string-method config unless overridden
if [ -z "$NI" ]; then
  NI=$(awk '$1=="set" && $2=="num_replicas" {print $3; exit}' "$STRING_METHOD_CONFIG")
fi
[ -n "$NI" ] || { echo "ERROR: could not read num_replicas from $STRING_METHOD_CONFIG." >&2; exit 1; }
[ -z "$NP" ] && NP=$NI

if [ $((NP % NI)) -ne 0 ]; then
  echo "ERROR: total ranks ($NP) is not a multiple of replicas ($NI)." >&2
  echo "       NAMD divides ranks evenly across replica partitions." >&2
  exit 1
fi

[ -d "$OUTDIR" ] || "$CFACV_BASEDIR/sbin/make_output_dirs.sh" "$OUTDIR" "$NI"

export CFACV_BASEDIR
echo "runstring.sh: $NI replicas, $NP ranks ($((NP / NI)) per replica)"
echo "runstring.sh: cfacv  $CFACV_BASEDIR"
echo "runstring.sh: namd3  $(command -v "$NAMD3")"

exec $MPIRUN -np "$NP" $MPIRUN_OPTS "$NAMD3" \
     +replicas "$NI" "$CONF" +stdout "$OUTDIR/%d/${CONF%.conf}.%d.log"
