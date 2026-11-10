#!/bin/bash
#SBATCH --job-name=cfacv-spr-sweep
#SBATCH --partition=def-sm
#SBATCH --account=abramsfreeprj
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Replica-synchronization sweep for the cfacv string method under NAMD 3.
#
# Holds the TOTAL number of MD timesteps fixed and varies steps_per_run, the
# number of MD steps between string updates.  Fixed total work means the only
# thing changing across runs is how many times the 24 replicas synchronize:
#
#     T(spr) = (S/spr) * t_sync  +  S * t_md
#
# so T plotted against 1/spr is a straight line whose slope is S*t_sync and
# whose intercept is S*t_md.  That separates barrier cost from MD cost
# without needing either to be measured directly.
#
# The initial string is reused from job 26741168 rather than regenerated, so
# every point starts from byte-identical coordinates and restraints.
set -euo pipefail

module purge
module load abramsGrp-namd/3.0.3/mpi-gcc

SRC=$SLURM_SUBMIT_DIR
CFACV=/ifs/groups/abramsGrp/pg565/cfacv
EXAMPLE=$CFACV/examples/alad-string-method
WORK=/scratch/$USER/$SLURM_JOB_ID
RESULTS=$SRC/sweep.$SLURM_JOB_ID

TOTAL_STEPS=20000
SPRS="1 2 5 10 20 50 100"

mkdir -p "$WORK" "$RESULTS"
CSV=$WORK/sweep.csv
echo "steps_per_run,num_runs,total_md_steps,wallclock_s,ms_per_iteration" > "$CSV"

stage_out() {
  cp -a "$CSV" "$RESULTS"/ 2>/dev/null || true
  for s in $SPRS; do
    [ -d "$WORK/spr_$s" ] || continue
    mkdir -p "$RESULTS/spr_$s"
    cp -a "$WORK/spr_$s"/output/0/*.history "$RESULTS/spr_$s"/ 2>/dev/null || true
    cp -a "$WORK/spr_$s"/output/0/job0.0.log "$RESULTS/spr_$s"/ 2>/dev/null || true
    cp -a "$WORK/spr_$s"/run.log "$RESULTS/spr_$s"/ 2>/dev/null || true
  done
}
trap stage_out EXIT

for spr in $SPRS; do
  nr=$(( TOTAL_STEPS / spr ))
  d=$WORK/spr_$spr
  mkdir -p "$d"
  cp -a "$EXAMPLE"/. "$d"/
  cp -a "$SRC"/initial_string/. "$d"/
  cd "$d"
  rm -rf output

  sed -i "s|^set steps_per_run .*|set steps_per_run $spr  ; # sweep|" alad_stringmethod.conf
  sed -i "s|^set num_runs .*|set num_runs $nr  ; # sweep: total MD steps held at $TOTAL_STEPS|" alad_stringmethod.conf

  echo "=== steps_per_run=$spr  num_runs=$nr  total=$TOTAL_STEPS"
  ./runstring.sh --cfacv-base-dir="$CFACV" --num-procs="$SLURM_NTASKS" > run.log 2>&1

  wc_s=$(grep -h "WallClock:" output/0/job0.0.log | tail -1 | awk "{print \$2}")
  if [ -z "$wc_s" ]; then
    echo "  WARNING: no WallClock for spr=$spr; see $d/run.log"
    echo "$spr,$nr,$TOTAL_STEPS,NA,NA" >> "$CSV"
    continue
  fi
  ms=$(awk -v t="$wc_s" -v n="$nr" "BEGIN{printf \"%.3f\", 1000*t/n}")
  echo "  WallClock ${wc_s}s   ${ms} ms/iteration"
  echo "$spr,$nr,$TOTAL_STEPS,$wc_s,$ms" >> "$CSV"
done

echo "=== sweep complete"
cat "$CSV"
