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
CFACV_BASEDIR=/home/cfa/research/cfacv
CHARMRUN=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/charmrun
NAMD2=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/namd2 
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
    CHARMRUN="${i#*=}"
    shift # past argument=value
    ;;
    -namd2=*|--namd2-command=*)
    NAMD2="${i#*=}"
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
for f in $CVINP $RESTRX $INITX $CHARMRUN $NAMD2 $STRING_METHOD_CONFIG; do
   if [ ! -f $f ]; then
      echo "ERROR: File $f not found."
      exit
   fi
done
for d in $CFACV_BASEDIR; do
   if [ ! -d $d ]; then
      echo "ERROR: Directory $d not found."
      exit
   fi
done

rootdir=`pwd`

N_SYSTEMS=`grep num_replicas $STRING_METHOD_CONFIG|awk '{print $3}'`         ; # number of actual MD systems NAMD/replica will handle
DUAL=`grep "SMPARAMS(dual)" $STRING_METHOD_CONFIG|awk '{print $3}'`
SYSTEMS_PER_IMAGE=1
if [ "$DUAL" -eq "1" ] ; then
  SYSTEMS_PER_IMAGE=2
fi

N_IMAGES=`echo "$N_SYSTEMS / $SYSTEMS_PER_IMAGE" | bc`
K=(40.0 1000.0)     ; # spring constants for restrained MD warmups stage-1 and stage-2
N_WARMUP_STAGES=${#K[@]}

N_RUNS=`echo "$N_WARMUP_STAGES * $N_SYSTEMS" | bc` ; # total number of MD simulations

echo "./mkinitimages will conduct $N_RUNS MD simulations"

NII=`echo "$N_IMAGES - 1" | bc`
NSI=`echo "$SYSTEMS_PER_IMAGE - 1"|bc`

echo "nproc $NP"
echo "nimg $N_IMAGES"
echo "syspi $SYSTEMS_PER_IMAGE"
echo "nsys $N_SYSTEMS"
echo "nwus $N_WARMUP_STAGES"

# initial string endpoints; string will be a line between these two points
# PHIPSI1=(-0.5 2.26849)
PHIPSI1=(1.8 -0.5)
PHIPSI2=(-0.25 -1.3)

DELPHI=`echo "scale=5; (${PHIPSI2[0]})-(${PHIPSI1[0]})"|bc -l`
DELPSI=`echo "scale=5; (${PHIPSI2[1]})-(${PHIPSI1[1]})"|bc -l`
DPHI=`echo "scale=5; $DELPHI / $NII"|bc -l`
DPSI=`echo "scale=5; $DELPSI / $NII"|bc -l`

# generate initial image systems
run=0
#for img in 0 ; do
for img in `seq 0 $NII`; do
   PHI=`echo "scale=5; (${PHIPSI1[0]}) + $img * $DPHI"|bc -l`
   PSI=`echo "scale=5; (${PHIPSI1[1]}) + $img * $DPSI"|bc -l` 
   for s in `seq 0 $NSI`; do
      system=`echo "$img + $s * $N_IMAGES"|bc`
      for i in `seq 1 $N_WARMUP_STAGES`; do
        ii=`echo "$i - 1"|bc`
        rfn=init_restr_${system}-${i}.inp
        if [ $i -eq $N_WARMUP_STAGES ] ; then
          rfn=restr.job0.${system}.inp
        fi
        cat restr_initX.inp | sed s/%K%/${K[$ii]}/g | sed s/%PHI%/$PHI/g | sed s/%PSI%/$PSI/g > $rfn
        cat initX.conf | \
            sed s/%IMG%/$img/g | \
            sed s/%REP%/$system/g | \
            sed s/%STAGE%/$i/g  | \
            sed s/%SEED%/$RANDOM/ | \
            sed s/%NSTAGES%/$N_WARMUP_STAGES/ > init_${system}-${i}.conf
        run=`echo "$run + 1"|bc`
        echo "Run $run img $img system $system stage $i in `pwd`"
        ${CHARMRUN} -n $NP ${NAMD2} init_${system}-${i}.conf > init_${system}-${i}.log
      done
   done
done

echo "Done."

