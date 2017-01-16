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
# (c) 2016 Cameron F Abrams, Drexel University
#
if [ ! -f cv.inp ]; then
  echo "ERROR: No cv.inp found."
  exit
fi
if [ ! -f restr_initX.inp ]; then
  echo "ERROR: No restr_initX.inp found."
  exit
fi
if [ ! -f initX.conf ]; then
  echo "ERROR: No initX.conf found."
  exit
fi

###################################################################
#### edit these directory/file names to conform to your system ####
CFACV_BASEDIR=/home/cfa/research/cfacv
CHARMRUN=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/charmrun
NAMD2=/home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/namd2 
###################################################################

NP=16 ; # number of processors to use for each MD simulation
#CWD=`pwd`
#NCV=`grep -v ^\# cv.inp | wc -l | awk '{print $1}'`
#CVL=`grep -v ^\# cv.inp | awk '{print $1}'`


NI=24 ; # number of images to use in string
NII=`echo "$NI - 1" | bc`

# initial string endpoints; string will be a line between these two points
PHIPSI1=(-0.5 2.26849)
PHIPSI2=(-0.25 -1.3)

DELPHI=`echo "scale=5; (${PHIPSI2[0]})-(${PHIPSI1[0]})"|bc -l`
DELPSI=`echo "scale=5; (${PHIPSI2[1]})-(${PHIPSI1[1]})"|bc -l`
DPHI=`echo "scale=5; $DELPHI / $NII"|bc -l`
DPSI=`echo "scale=5; $DELPSI / $NII"|bc -l`

#echo "DELPHI $DELPHI DELPSI $DELPSI DPHI $DPHI DPSI $DPSI NI $NI"

# set up and run sequential low-k (40) and high-k (100) restrained MD simulations
# to generate initial image systems
for img in `seq 0 $NII`; do
   PHI=`echo "scale=5; (${PHIPSI1[0]}) + $img * $DPHI"|bc -l`
   PSI=`echo "scale=5; (${PHIPSI1[1]}) + $img * $DPSI"|bc -l` 
   i=1
   for k in 40.0 1000.0; do
      cat restr_initX.inp | sed s/%K%/$k/g | sed s/%PHI%/$PHI/g | sed s/%PSI%/$PSI/g > restr_init_${img}-${i}.inp
      cat initX.conf | sed s/%IMG%/$img/g | sed s/%STEP%/$i/g > init_${img}-${i}.conf
      echo "Running img $img step $i"
      ${CHARMRUN} -n $NP ${NAMD2} init_${img}-${i}.conf > init_${img}-${i}.log
      i=`echo "$i + 1"|bc`
   done
done

puts "Done."

