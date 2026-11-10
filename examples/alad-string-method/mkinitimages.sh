#!/bin/bash
#
# String method in collective variables via NAMD/replica
#
# Creation of initial images
#
# This script sets up initial images of alanine dipeptide
# in (phi,psi) space for a string method calculation.
#
# Required files:
#   - cv.inp used by cfacv.tcl in the TAMD sweeps
#   - restr_initX.tmp - template restr.inp file for restrained MD via cfacv
#   - initX.conf - template NAMD configuration file for the restrained MD
#
# (c) 2016-2018 Cameron F Abrams, Drexel University
# cfa22@drexel.edu
#
# default input file names
CVINP="cv.inp"
RESTRX="restr_initX.inp"
INITX="initX.conf"
# Reworked 2026-11-10: resolve cfacv and NAMD from the environment rather
# than from hardcoded paths on the original author's workstation.  These are
# plain restrained-MD warmups, one system at a time and no replicas, so the
# single-node build is the right one and charmrun is not involved.
#
# With namd3 on PATH:
#
#     ./mkinitimages.sh --num-procs=8
#
CFACV_BASEDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
NAMD3=namd3
NAMD_CMD=""      # full launcher; default is "$NAMD3 +p$NP"
NP=16
STRING_METHOD_CONFIG="alad_stringmethod.conf"

for i in "$@"
do
case $i in
    -cvi=*|--cv-input-file=*)
    CVINP="${i#*=}"
    shift # past argument=value
    ;;
    -rstx=*|--restraint-template-file=*)
    RESTRX="${i#*=}"
    shift # past argument=value
    ;;
    -inx=*|--init-template-file=*)
    INITX="${i#*=}"
    shift # past argument=value
    ;;
    -cfacvbd=*|--cfacv-base-dir=*)
    CFACV_BASEDIR="${i#*=}"
    shift # past argument=value
    ;;
    -charmrun=*|--charmrun-command=*)
    echo "mkinitimages.sh: --charmrun is obsolete and ignored; the warmups are"
    echo "                 single-system runs launched as '$NAMD3 +p<NP>'."
    shift # past argument=value
    ;;
    -namd3=*|--namd3-command=*|-namd2=*|--namd2-command=*)
    NAMD3="${i#*=}"
    shift # past argument=value
    ;;
    --namd-cmd=*)
    NAMD_CMD="${i#*=}"
    shift # past argument=value
    ;;
    -np=*|--num-procs=*)
    NP="${i#*=}"
    shift # past argument=value
    ;;
    -i=*|--string-method-config=*)
    STRING_METHOD_CONFIG="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
          # unknown option
    ;;
esac
done
for f in $CVINP $RESTRX $INITX $STRING_METHOD_CONFIG; do
   if [ ! -f $f ]; then
      echo "ERROR: File $f not found."
      exit 1
   fi
done
if [ ! -d "$CFACV_BASEDIR/tcl" ]; then
   echo "ERROR: $CFACV_BASEDIR is not a cfacv tree."
   exit 1
fi
if [ -z "$NAMD_CMD" ]; then
   NAMD_CMD="$NAMD3 +p$NP"
fi
if ! command -v ${NAMD_CMD%% *} >/dev/null; then
   echo "ERROR: '${NAMD_CMD%% *}' not on PATH."
   echo "       module load abramsGrp-namd/3.0.3/multicore-gcc"
   exit 1
fi
echo "mkinitimages.sh: cfacv $CFACV_BASEDIR"
echo "mkinitimages.sh: namd  $NAMD_CMD"
export CFACV_BASEDIR

rootdir=`pwd`

N_SYSTEMS=`grep num_replicas $STRING_METHOD_CONFIG|awk '{print $3}'`         ; # number of actual MD systems NAMD/replica will handle
DUAL=`grep "SMPARAMS(dual)" $STRING_METHOD_CONFIG|awk '{print $3}'`
SYSTEMS_PER_IMAGE=1
if [ "$DUAL" -eq "1" ] ; then
  SYSTEMS_PER_IMAGE=2
fi

N_IMAGES=$(( N_SYSTEMS / SYSTEMS_PER_IMAGE ))
K=(40.0 1000.0)     ; # spring constants for restrained MD warmups stage-1 and stage-2
N_WARMUP_STAGES=${#K[@]}

N_RUNS=$(( N_WARMUP_STAGES * N_SYSTEMS )) ; # total number of MD simulations

echo "./mkinitimages will conduct $N_RUNS MD simulations"

NII=$(( N_IMAGES - 1 ))
if [ "$NII" -lt 1 ]; then
   echo "ERROR: need at least 2 images; num_replicas=$N_SYSTEMS with"
   echo "       SMPARAMS(dual)=$DUAL gives $N_IMAGES image(s)."
   exit 1
fi
NSI=$(( SYSTEMS_PER_IMAGE - 1 ))

echo "nproc $NP"
echo "nimg $N_IMAGES"
echo "syspi $SYSTEMS_PER_IMAGE"
echo "nsys $N_SYSTEMS"
echo "nwus $N_WARMUP_STAGES"

# initial string endpoints; string will be a line between these two points
# PHIPSI1=(-0.5 2.26849)
PHIPSI1=(1.8 -0.5)
PHIPSI2=(-0.25 -1.3)

# Floating-point arithmetic via awk.  This used to shell out to bc, which is
# not installed on every host (it is absent from the group workstation), and
# its absence was silent: bc printed "command not found" to stderr, every
# arithmetic result came back empty, and the script still reported "Done."
fcalc() { awk "BEGIN{printf \"%.5f\", $1}"; }

DELPHI=$(fcalc "(${PHIPSI2[0]})-(${PHIPSI1[0]})")
DELPSI=$(fcalc "(${PHIPSI2[1]})-(${PHIPSI1[1]})")
DPHI=$(fcalc "$DELPHI / $NII")
DPSI=$(fcalc "$DELPSI / $NII")

# generate initial image systems
run=0
#for img in 0 ; do
for img in `seq 0 $NII`; do
   PHI=$(fcalc "(${PHIPSI1[0]}) + $img * $DPHI")
   PSI=$(fcalc "(${PHIPSI1[1]}) + $img * $DPSI")
   for s in `seq 0 $NSI`; do
      system=$(( img + s * N_IMAGES ))
      for i in `seq 1 $N_WARMUP_STAGES`; do
        ii=$(( i - 1 ))
        rfn=init_restr_${system}-${i}.inp
        if [ $i -eq $N_WARMUP_STAGES ] ; then
          rfn=restr.job0.${system}.inp
        fi
        cat restr_initX.inp | sed s/%K%/${K[$ii]}/g | sed s/%PHI%/$PHI/g | sed s/%PSI%/$PSI/g > $rfn
        cat initX.conf | \
            sed s,%CFACVBD%,$CFACV_BASEDIR,g | \
            sed s/%IMG%/$img/g | \
            sed s/%REP%/$system/g | \
            sed s/%STAGE%/$i/g  | \
            sed s/%SEED%/$RANDOM/ | \
            sed s/%NSTAGES%/$N_WARMUP_STAGES/ > init_${system}-${i}.conf
        run=$(( run + 1 ))
        echo "Run $run img $img system $system stage $i in `pwd`"
        ${NAMD_CMD} init_${system}-${i}.conf > init_${system}-${i}.log
      done
   done
done

echo "Done."

