#CFACV -- a collective variables implementation for NAMD

Cameron F Abrams, Drexel University, Philadelphia, Pennsylvania  
cfa22@drexel.edu  
215 895 2231  
(c) 2016   

###INTRODUCTION

CFACV is a simple collective variables implementation for NAMD that
uses 'tclforces' and implements Temperature-Accelerated Molecular
Dynamics (TAMD).  Collective variables are functions of atomic
Cartesian coordinates, and CFACV currently allows for CV's to be
centers of mass of groupings of atoms and any lengths, angles, and
dihedrals constructed thereof.  These CV's can be harmonically
restrained or driven with Temperature-Accelerated MD.

In this friendly release, there are three subdirectories, summarized below:

1. src -- C files and headers, plus a makefile
2. tcl -- tcl scripts
3. example -- an example TAMD run


###INSTRUCTIONS FOR RUNNING THE HIV-1 MA EXAMPLE

1. Compile shared object libraries cfacv.so and libgenericdataspace.so

   From within the main directory:
    * mkdir lib
    * cd src
    * make cfacv.so
    * make libgenericdataspace.so 

     Note that you must have swig installed
 
2. Run the example (This is a system with HIV-1 MA (PDB 2hiw))

   * cd into examples/hiv1-ma
   * examine mk_tPDB.tcl; you can edit parts of this
     Currently, this is set to block the protein into 5 subdomains
   * issue "vmd -dispdev text -e mk_tPDB.tcl -args -p my_system.pdb"
     This will create "label.pdb"; "label.pdb.BAK" is a copy I created.
   * issue a command to launch namd2 with the config file "example.conf"
     For example,
         charmrun ++local +p4 /usr/local/bin/namd2 example.conf +idlepoll
   * watch it run -- the log file has lines that begin with "CFACV/C) that
     report data on fictitious variables "Z" and congruent collective variables
     "Th" (for "theta")

The other example is alanine dipeptide in TIP3P water, with the phi and psi
angles as CV's, driven by TAMD.

###INSTRUCTIONS FOR INSTALLING CFACV

1. copy the cfacv/ directory wholesale under your home directory

2. in any run directory where you want to use cfacv, make sure the
   NAMD configuration file references the right path for cfacv_tclforces.tcl
   as being ${HOME}/cfacv/tcl/cfacv_tclforces.tcl

###MORE DETAILS

CFACV interfaces with NAMD via a tclforces script called
'cfacv_tclforces.tcl' that initializes a workspace and then defines
the all-important 'calcforces' procedure.  Procedures called from this
script are defined in cfacv.tcl and much of the implementation is
housed in a cfacv.c.  Swig is used to interface cfacv.tcl and cfacv.c
to create a shared object library 'cfacv.so' loaded by
cfacv_tclforces.  (Note: you need swig installed to build the
shared-object library.)

The configuration file 'example.conf' contains a tclforces stanza
for a TAMD simulation.  In addition to defining the tclforces script, the TAMD
implementations need the names of three accessory files:

* labelPDB: this is a PDB congruent to the simulation system in which
 the beta field of each atom designates its membership in a "group". (0
 designates the "null" group.)  You have to create this file yourself.
 The mk_tPDB.tcl script illustrates how one can use automatic clustering
 of residues into subdomains to generate subdomain groups.  Groups
 can have a minimum of one atom.

* cvINP: this is a custom input file that identifies the types of
 collective variables and what groups they involve.  The example
 cv.inp in this directory defines the cartesian coordinates of the
 center of mass of three groups as six collective variables.  There are
 other types of CV's; consult the source.

 NOTE:  the group ID explicitly listed in the beta field of the labelPDB
 is an index that begins at 1.  "1" in this context means "group 0", the first
 group.  Groups listed in the cvINP file indexed beginning at 0.  So, all atoms
 with a "1" in their beta field in the labelPDB correspond to "group 0" in the
 cvINP file; those with "2" in their beta fields are "group 1", etc.

* restrINP: this is a custom input file that defines how the CV's are
 restrained and/or driven by tclforces.  It has a matrix definition of
 the restrained variables (which are generally linear combinations of CV's),
 and each row also lists TAMD parameters in a self-explanatory way.
 NOTE:  this input file is optional; if it is not listed, you must
 have a line in the NAMD configuration file that looks something like

> set restrPARAMS    {{k 100} {TAMDkT 6.0} {TAMDgamma 100.0} {TAMDdt 0.002}}

 In either case, time is measured in *ps* for all TAMD variables.

To understand the implementation, I recommend tracing the code
beginning with cfacv_tclforces.tcl.  Issue 'make cfacv.so' to build the
shared object library.

If you are only going to accelerate Cartesian variables, this is good
enough.  If you are going to accelerate dihedral angles, you have to
indicate that the restraint function in periodic with the additional
TAMD parameter {rf PERIODIC}.

Please contact me with any questions or problems.  This is not a professional
software release, so I fully anticipate there are bugs or other confusing
issues.

