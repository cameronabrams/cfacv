#!/bin/bash
#
# String method in collective variables via NAMD/replica
#
# Measurement of squared displacement of a string from a reference
#
# (c) 2016 Cameron F Abrams, Drexel University
#
# default values
cv_inp=cv.inp
sm_history=output/0/alad_sm.job0.0.history
CFACV_BASEDIR=/home/cfa/research/cfacv

avg_last=1  # number of iterations to average over to generate anchor points
NI=24  # number of images
DIHED="-dihed" # set to "-dihed" if CV's are dihedral angles

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -c|--cv_inp)
    cv_inp="$2"
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
    -ni|--number_of_images)
    NI="$2"
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

for f in $cv_inp $string_method_history ; do
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

NCV=`grep -v ^\# $cv_inp | wc -l | awk '{print $1}'`
NCVI=`echo "$NCV - 1" | bc`
CVL=`grep -v ^\# $cv_inp | awk '{print $1}'`
NII=`echo "$NI - 1" | bc`

NLH=`grep -w reparam $sm_history | wc -l | awk '{print $1}'`
NITER=`echo "scale=0; $NLH / $NI" | bc -l`
echo "# String euclidean distance from reference string, averaged over images"
echo "# Extracted from history file $sm_history which contains info on $NITER iterations."
echo "# Extracting CV values from last $avg_last iterations and averaging to produce reference string"
echo "# Data follows: [iteration] [RMSD]"
#echo "Initializing Z..."
for j in `seq 0 $NII`; do
   for i in `seq 0 $NCVI`; do
     idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
     Z[$idx]=0
#     echo "  INIT: img $j cv-cmp $i idx $idx Z[$idx] ${Z[$idx]}"
   done
done
#echo "Visiting iterations..."
for nn in `seq 1 $avg_last`; do
#   echo "Iteration $nn..."
   NLL=`echo "scale=0; ${NI} * $nn" | bc -l`
   for i in `seq 0 $NCVI`; do
#     echo "Iteration $nn CV component $i..."
     ii=`echo "$i + 1" | bc`
     THISZ=(`grep -w reparam $sm_history | tail -${NLL} | head -${NI} | awk -F'z ' '{print $2}' | cut -d' ' -f $ii`)
#     echo -n "   -> per-image Z "
     for j in `seq 0 $NII`; do
#       echo -n "${THISZ[$j]} "
       idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
       Z[$idx]=`echo "scale=8; ${Z[$idx]} + ${THISZ[$j]}" | bc -l`
     done
#     echo "."
   done
done
#echo "Dividing tallies..."
for j in `seq 0 $NII`; do
   for i in `seq 0 $NCVI`; do
     idx=`echo "scale=0; $j * $NCV + $i"|bc -l`
     Z0[$idx]=`echo "scale=8; ${Z[$idx]} / $avg_last" | bc -l`
   done
done
#echo "Completed data averaging."

# get Z for each frame and compute squared euclidean distance from Z0
for nn in `seq 1 $NITER`; do
#   echo "Iteration $nn..."
   NLL=`echo "scale=0; ${NI} * $nn" | bc -l`
   for i in `seq 0 $NCVI`; do
#     echo "Iteration $nn CV component $i..."
     ii=`echo "$i + 1" | bc`
     THISZ=(`grep -w reparam $sm_history | head -${NLL} | tail -${NI} | awk -F'z ' '{print $2}' | cut -d' ' -f $ii`)
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
   echo "$nn $RMSD"
done

#echo "Done."
