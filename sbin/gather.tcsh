#!/bin/tcsh
#
# Single-sweep free-energy reconstruction
#
# Mean-force calculations on centers extracted from tamd sweeps
# via cfacv_harvest_centers.tcl
#
# This script optionally serially runs all restrained force-measuring
# MD simulations for each centerX.coor and then processes each
# resulting log file to extract the mean forces.
#
# Required files:
#   - cv.inp used by cfacv.tcl in the TAMD sweeps
#   - mesh.dat created by the cfacv_harvest_centers.tcl
#   - restrfmX.tmp - template restr.inp file for restrained MD via cfacv
#     this script replaces CV1 with the value of 1st component of the CV restraint, etc.
#   - fmX.conf - template NAMD configuration file for the restrained MD
#
# Generates output forces.dat -- formatted for use by Reconstruct.c
#
# (c) 2016 Cameron F Abrams, Drexel University
#
if ( ! -f mesh.dat ) then
  echo "ERROR: No mesh.dat found.  Did you harvest points yet?"
  exit
endif
if ( ! -f cv.inp ) then
  echo "ERROR: No cv.inp found."
  exit
endif
if ( ! -f restrfmX.inp ) then
  echo "ERROR: No restrfmX.inp found."
  exit
endif
if ( ! -f fmX.conf ) then
  echo "ERROR: No fmX.conf found."
  exit
endif
setenv REGENERATE 1 # set to 1 if you want to force regeneration via MD
setenv DIHED -dihed # set to 1 if CV's are dihedral angles on domain [-pi:pi]
setenv CFACV_BASEDIR /home/cfa/research/cfacv
setenv CHARMRUN /home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/charmrun 
setenv NP 6  # number of processors to use for each MD simulation
setenv CWD `pwd`
setenv NCV `grep -v ^\# cv.inp | wc -l | awk '{print $1}'`
setenv CVL `grep -v ^\# cv.inp | awk '{print $1}'`
setenv NPOINTS `grep -v ^\# mesh.dat | wc -l | awk '{print $1}'`
setenv NAMD2 /home/cfa/namd/NAMD_2.11_Source/Linux-x86_64-g++/namd2
cat > forces.dat << EOF
# forces.dat
# $CWD
#
$NCV $NPOINTS
$CVL

EOF
set m = 1
foreach d (center*.coor)
  set cv = (`head -$m mesh.dat | tail -1`)
  set n = `echo $d | sed s/center//`
  set n = `echo $n | sed s/.coor//`
  echo "$n"
  if ( ! -f fm${n}.log || $REGENERATE ) then
    set i = 1
    cp restrfmX.inp tmp
    while ( $i <= $NCV) 
      cat tmp | sed s/%CV${i}%/$cv[$i]/g > tmp2
      mv tmp2 tmp
      @ i++
    end
    mv tmp restrfm${n}.inp
    cat fmX.conf | sed s/%X%/$n/g > fm${n}.conf
    ${CHARMRUN} -n $NP ${NAMD2} fm${n}.conf > fm${n}.log
  endif
  ${CFACV_BASEDIR}/bin/forcesFromLog -f fm${n}.log -dim $NCV -ofn fra${n}.dat $DIHED >> forces.dat
  @ m++
end

tail -$m forces.dat > forces_gnuplot.dat

