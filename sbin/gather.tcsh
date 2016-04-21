#!/bin/tcsh
#
# This script gathers forces from the log files in each
# Frame${N} directory.
#
# (c) 2010 Cameron F Abrams, Drexel University
#
set last = "NONE"
if ( -f forces.dat)  then
  set last = `tail -1 forces.dat | awk '{print $1}'`
else
  echo "# CV    CV forces"
endif
foreach d (Frame*)
  cd $d
  set n = `echo $d | sed s/Frame//`
  if ( -f sample.log ) then
      set dat=`${HOME}/bin/forcesFromLog -f sample.log`
      echo "$n $dat"
  endif
  cd ..
end
