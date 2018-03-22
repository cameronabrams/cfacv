#!/bin/bash
#
# String method in collective variables via NAMD/replica
#
# Measurement of string RMSD
#
# Required files:
#   - cv.inp used by cfacv.tcl in the TAMD sweeps
#   - restr_initX.tmp - template restr.inp file for restrained MD via cfacv
#   - initX.conf - template NAMD configuration file for the restrained MD
#   - stringmethod.conf - string method configuration file used in the SMCV simulation
#
# (c) 2016-2018 Cameron F Abrams, Drexel University
#
# default values
cv_inp=cv.inp
string_method_history=output/0/alad_sm.job0.0.history
STRING_METHOD_CONFIG=alad_stringmethod.conf
CFACV_BASEDIR=/home/cfa/research/cfacv

avg_last=1  # number of iterations to average over to generate anchor points
DIHED="-dihed" # set to "-dihed" if CV's are dihedral angles

while [[ $# -gt 1 ]]
do
key="$1"
case $key in
    -c|--cv_inp)
    cv_inp="$2"
    shift # past argument
    ;;
    -smc|--string_method_config)
    STRING_METHOD_CONFIG="$2"
    shift # past argument
    ;;
    -smh|--string_method_history)
    string_method_history="$2"
    shift # past argument
    ;;
    -nA|--avg_last)
    avg_last="$2"
    shift # past argument
    ;;
    -cbd|--cfacv_basedir)
    CFACV_BASEDIR="$2"
    shift
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

for f in $cv_inp $STRING_METHOD_CONFIG $string_method_history ; do
  if [ ! -f $f ]; then
    echo "ERROR: file $f not found."
    exit
  fi
done

for d in $CFACV_BASEDIR; do
  if [ ! -d $d ]; then
    echo "ERROR: directory $d not found."
    exit
  fi
done

N_SYSTEMS=`grep num_replicas $STRING_METHOD_CONFIG|awk '{print $3}'`
DUAL=`grep "SMPARAMS(dual)" $STRING_METHOD_CONFIG|awk '{print $3}'`
SYSTEMS_PER_IMAGE=1
if [ "$DUAL" -eq "1" ] ; then
  SYSTEMS_PER_IMAGE=2
fi

NI=`echo "$N_SYSTEMS / $SYSTEMS_PER_IMAGE" | bc`
NCV=`grep -v ^\# $cv_inp | wc -l | awk '{print $1}'`
NCVI=`echo "$NCV - 1" | bc`
CVL=`grep -v ^\# $cv_inp | awk '{print $1}'`
NII=`echo "$NI - 1" | bc`

NLH=`grep -w reparam $string_method_history | wc -l | awk '{print $1}'`
NITER=`echo "scale=0; $NLH / $N_SYSTEMS" | bc -l`
echo "History file $string_method_history contains info on $NITER iterations over $N_SYSTEMS systems and $NI images."

echo "Extracting CV values from last $avg_last iterations"
echo "Initializing Z..."
for j in `seq 0 $NII`; do
   for i in `seq 0 $NCVI`; do
     idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
     Z[$idx]=0
#     echo "  INIT: img $j cv-cmp $i idx $idx Z[$idx] ${Z[$idx]}"
   done
done
echo "Visiting iterations..."
for nn in `seq 1 $avg_last`; do
   echo "Iteration $nn..."
   NLL=`echo "scale=0; ${N_SYSTEMS} * $nn" | bc -l`
   for i in `seq 0 $NCVI`; do
     echo "Iteration $nn CV component $i..."
     ii=`echo "$i + 1" | bc`
     THISZ=(`grep -w reparam $string_method_history | tail -${NLL} | head -${NI} | awk -F'z ' '{print $2}' | cut -d' ' -f $ii`)
     echo -n "   -> per-image Z "
     for j in `seq 0 $NII`; do
       echo -n "${THISZ[$j]} "
       idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
       Z[$idx]=`echo "scale=8; ${Z[$idx]} + ${THISZ[$j]}" | bc -l`
     done
     echo "[END]"
   done
done

echo "Dividing tallies..."
for j in `seq 0 $NII`; do
   for i in `seq 0 $NCVI`; do
     idx=`echo "scale=0; $j * $NCV + $i"|bc -l`
     Z0[$idx]=`echo "scale=8; ${Z[$idx]} / $avg_last" | bc -l`
   done
done
echo "Completed data averaging."

# get Z for each frame and compute squared euclidean distance from Z0
if [ -f rmsd.dat ] ; then
  echo "Removing rmsd.dat"
fi
echo "# iter RMSD" > rmsd.dat

for nn in `seq 1 $NITER`; do
#   echo "Iteration $nn..."
   NLL=`echo "scale=0; ${N_SYSTEMS} * $nn" | bc -l`
   for i in `seq 0 $NCVI`; do
#     echo "Iteration $nn CV component $i..."
     ii=`echo "$i + 1" | bc`
     THISZ=(`grep -w reparam $string_method_history | head -${NLL} | tail -${NI} | awk -F'z ' '{print $2}' | cut -d' ' -f $ii`)
#     echo -n "   -> per-image Z "
     for j in `seq 0 $NII`; do
#       echo -n "${THISZ[$j]} "
       idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
       Z[$idx]=`echo "scale=8; ${Z0[$idx]} - ${THISZ[$j]}" | bc -l`
       Z[$idx]=`echo "scale=8; ${Z[$idx]} * ${Z[$idx]}"   | bc -l`
     done
#     echo "."
   done
   RMSD=0
   for j in `seq 0 $NII`; do
     D[$j]=0
     for i in `seq 0 $NCVI`; do
       idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
       D[$j]=`echo "scale=8; ${D[$j]} + ${Z[$idx]}" | bc -l`
     done
     RMSD=`echo "scale=8; $RMSD + ${D[$j]}" | bc -l`
   done
   RMSD=`echo "scale=8; sqrt($RMSD / $NI)" | bc -l`
   echo "$nn $RMSD" >> rmsd.dat
   echo $nn
done

echo "Done.  Created rmsd.dat."
