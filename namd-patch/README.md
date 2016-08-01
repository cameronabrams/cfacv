#The Position-Velocity-Rewind patch for NAMD 2.11

Cameron F Abrams, Drexel University, Philadelphia, Pennsylvania  
cfa22@drexel.edu  
215 895 2231  
(c) 2016   


PVRW is a patch to the NAMD 2.11 source code that implements
the capability, on demand during a running MD simulation,
rewind all atomic position by one time step AND rewind
all atomic velocities by one time step and negate them.


To use PVRW, copy the patch file to the NAMD_2.11_Source/ directory
and then patch the src/

patch -p0 < cfa_pvrw_namd211.patch

Now compile NAMD as you would normally.

Next, you must issue a command of the form

setboundaryflag $step 0 

every time step inside the calc_forces routine of tcl forces.  If a condition
is met where you would like to execute a PVRW, you instead issue

setboundaryflag $step 1

where in both cases, "step" is the current timestep.

Note that this command has to be executed at every time step with either
 a "0" or "1" as the second argument since it sets a broadcast that
the sequencer expects.  The code will wait patiently for this broadcast 
if you don't send it!

