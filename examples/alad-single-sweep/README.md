#Single-sweep reconstruction for computing free energies

Cameron F Abrams, Drexel University, Philadelphia, Pennsylvania  
cfa22@drexel.edu  
215 895 2231  
(c) 2016   

This example uses NAMD and the CFACV collective variables module to implement
the method of single-sweep reconstruction to compute the free energy
in the space of (phi,psi) of alanine dipeptide, replicating
the results in Maragliano and Vanden-Eijnden, JCP 128, 184110 (2008).


###INSTRUCTIONS

1. Make sure everything in $CFACV_BASEDIR/src is compiled; "make all"

2. Conduct a TAMD "sweep" using NAMD; see tamd.conf

3. Harvest centers using $CFACV_BASEDIR/tcl/cfacv_harvest_centers.tcl

vmdtext -e ../../tcl/cfacv_harvest_centers.tcl -args -psf alad_wb.psf -dcd go_tamd.dcd -cv cv.inp -restr restr.inp -stip 0.373

4. Use the gather.tcsh script in sbin/ to setup, run, and analyze 
the restrained MD runs at each center

../../sbin/gather.tcsh > & gather.log &

5. Use Reconstruct.c to first search for the optimum value of sigma

../../bin/Reconstruct -sigscan 0.3,0.8,0.01  > & sigscan.out &

6.  Use Reconstruct.c again, at the optimum sigma value, to generate the final reconstruction and map:

../../bin/Reconstruct -sig 0.64 -mapres 0.0087

7. Plot using gnuplot script plotmap.gp

Note that this uses the file ref_md.dist, which is a probability distribution and free energy
extracted from a direct MD simulation.

