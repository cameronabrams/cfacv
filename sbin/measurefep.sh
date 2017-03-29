#!/bin/bash
#
# String method in collective variables via NAMD/replica
#
# Measurement of FEP along string
#
# This script sets up initial images of alanine dipeptide
# in (phi,psi) space for a string method calculation.
#
# Required files:
#   - cv.inp used by cfacv.tcl in the TAMD sweeps
#   - restr_initX.tmp - template restr.inp file for restrained MD via cfacv
#   - initX.conf - template NAMD configuration file for the restrained MD
#
# (c) 2016 Cameron F Abrams, Drexel University
#
# default values
cv_inp=cv.inp
restr_inp_template=restr_fepX.inp
namd_config_template=fepX.conf
sm_history=output/0/alad_sm.job0.0.history
CFACV_BASEDIR=/home/cfa/research/cfacv
CHARMRUN=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/charmrun
NAMD2=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/namd2

avg_last=1  # number of iterations to average over to generate anchor points
NI=24  # number of images
NP=16  # number of processors to use for each MD simulation
k=100.0  # spring constant in restrained MD
DIHED="-dihed" # set to "-dihed" if CV's are dihedral angles

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -c|--cv_inp)
    cv_inp="$2"
    shift # past argument
    ;;
    -rX|--restr_inp_template)
    restr_inp_template="$2"
    shift # past argument
    ;;
    -nX|--namd_config_template)
    namd_config_template="$2"
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
    -cr|--charmrun)
    CHARMRUN="$2"
    shift
    ;;
    -namd2)
    NAMD2="$2"
    shift
    ;;
    -np|--number_of_processors)
    NP="$2"
    shift
    ;;
    -ni|--number_of_images)
    NI="$2"
    shift
    ;;
    -k|--spring_constant)
    k="$2"
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

for f in $cv_inp $restr_inp_template $namd_config_template $string_method_history $CHARMRUN $NAMD2; do
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

image_forces_dat=image_forces.dat

NLH=`grep -w reparam $sm_history | wc -l | awk '{print $1}'`
NITER=`echo "scale=0; $NLH / $NI" | bc -l`
echo "History file $sm_history contains info on $NITER iterations."
echo "Extracting CV values from last $avg_last iterations"
echo "Initializing Z..."
for j in `seq 0 $NII`; do
   for i in `seq 0 $NCVI`; do
     idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
     Z[$idx]=0
     echo "  INIT: img $j cv-cmp $i idx $idx Z[$idx] ${Z[$idx]}"
   done
done
echo "Visiting iterations..."
for nn in `seq 1 $avg_last`; do
   echo "Iteration $nn..."
   NLL=`echo "scale=0; ${NI} * $nn" | bc -l`
   for i in `seq 0 $NCVI`; do
     echo "Iteration $nn CV component $i..."
     ii=`echo "$i + 1" | bc`
     THISZ=(`grep -w reparam $sm_history | tail -${NLL} | head -${NI} | awk -F'z ' '{print $2}' | cut -d' ' -f $ii`)
     echo -n "   -> per-image Z "
     for j in `seq 0 $NII`; do
       echo -n "${THISZ[$j]} "
       idx=`echo "scale=0; $j * $NCV + $i" | bc -l`
       Z[$idx]=`echo "scale=8; ${Z[$idx]} + ${THISZ[$j]}" | bc -l`
     done
     echo "."
   done
done
echo "Dividing tallies..."
for j in `seq 0 $NII`; do
   for i in `seq 0 $NCVI`; do
     idx=`echo "scale=0; $j * $NCV + $i"|bc -l`
     Z[$idx]=`echo "scale=8; ${Z[$idx]} / $avg_last" | bc -l`
   done
done
echo "Completed data averaging."

cat > $image_forces_dat << EOF
# image mean-forces computed using restrained MD with k = $k
# MD engine: $NAMD2
# ---z----   -----f-----
EOF
for img in `seq 0 $NII`; do
  cat $namd_config_template | sed s/%IMG%/$img/g > fep_${img}.conf
  rfn=`grep -w restrINP fep_${img}.conf|awk '{print $2}'`
  cat $restr_inp_template | sed s/%K%/$k/g > $rfn 
  echo -n "Image $img : Z "
  for i in `seq 0 $NCVI`; do
     idx=`echo "scale=0; $img * $NI + $i"|bc -l`
     echo -n "${Z[$idx]} "
     cp $rfn tmp; cat tmp | sed s/%CV${i}%/${Z[$idx]}/ > $rfn; rm tmp
  done
  echo " ...Running..."
  ${CHARMRUN} -n $NP ${NAMD2} fep_${img}.conf > fep_${img}.log
  ${CFACV_BASEDIR}/bin/ForcesFromLog -k $k -f fep_${img}.log -dim $NCV -ofn fep_fra_${img}.dat $DIHED >> $image_forces_dat
done

${CFACV_BASEDIR}/bin/FEPFromForces -i $image_forces_dat -dim $NCV -o fep.dat
echo "fep.dat created."

