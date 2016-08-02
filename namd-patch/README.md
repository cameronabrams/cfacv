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

Now, issue the ./config script with the following additional flags

--cxx-opts "-DCFA_PVRW" --cc-opts "-DCFA_PVRW"

Now, compile NAMD as you would normally.  The resulting
namd2 executable now has the PVRW capability enabled.

To use PVRW, it is assumed you are using tclforces to
set the PVRW flag.  That is, inside your calcforces
procedure, you are using the coordinates to
decide whether or not to execute a PVRW.  The way
the PVRW is implemented, you must set the pvrw
flag every time the calcforces procedure is called; that is, every time step.  This is because it is a "broadcast" (but I am sure I could figure out a way to do this
differently!).

As an example, consider this calcforces procedure:

  proc calcforces {} {
    set step [getstep]
    if {$step && ![expr $step % 100]} {
      puts "PVRW) Test: setpvrwflag $step 1"
      setpvrwflag $step 1
     } else {
      setpvrwflag $step 0
     }
  }

This procedure sets the PVRW flag to 1
if the timestep is a multiple of 100;
otherwise it sets the PVRW flag to 0.  Note that
the pvrw flag is explicitly set every time step!

The ala2-test/ directory contains all input files necessary
to run a simple test of NAMD 2.11 with PVRW enabled.  Two
representative logs are provided as well.

