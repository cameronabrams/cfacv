# CFACV -- a collective variables implementation for NAMD

Cameron F Abrams, Drexel University, Philadelphia, Pennsylvania  
cfa22@drexel.edu  
215 895 2231  
(c) 2016-2018   

CFACV is citable and archived at Zenodo: [![DOI](https://zenodo.org/badge/20736/cameronabrams/cfacv.svg)](https://zenodo.org/badge/latestdoi/20736/cameronabrams/cfacv)  

Please also cite "Cameron F. Abrams and Eric Vanden-Eijnden, "Large-scale conformational sampling of proteins using temperature-accelerated molecular dynamics," Proc. Natl. Acad. Sci. USA 107 4961-4966 (2010)"  

### INTRODUCTION

CFACV is a simple collective variables implementation for NAMD that
uses `tclforces` and implements Temperature-Accelerated Molecular
Dynamics (TAMD).  Collective variables are functions of atomic
Cartesian coordinates, and CFACV currently allows for CV's to be
centers of mass of groupings of atoms and any lengths, angles, and
dihedrals constructed thereof.  These CV's can be harmonically
restrained or driven with Temperature-Accelerated MD.

In this release (v1.0), there are three subdirectories, summarized below:

1. src -- C files and headers, plus a makefile
2. tcl -- tcl scripts
3. examples -- three example TAMD simulation inputs


### INSTRUCTIONS FOR RUNNING THE HIV-1 MA EXAMPLE

1. Compile shared object libraries cfacv.so and libgenericdataspace.so

   From within the main directory:
    * mkdir lib (you only have to do this once)
    * cd src
    * make cfacv.so
    * make genericdataspace.so

     Note that you must have swig installed.
 
2. Run the example (This is a system with HIV-1 MA (PDB 1hiw))

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

The example in the directory alad/ is CHARM22 alanine dipeptide (residue ALAD) in TIP3P water, with the phi and psi
angles as CV's, driven by TAMD.  The directory ala2-fixed-phi-psi/ is an example in which 
cfacv is used to simply restrain at a particular (phi,psi).

### INSTRUCTIONS FOR INSTALLING CFACV

Prerequisites:

1. swig and tcl headers and libraries

2. Gnu Scientific Library

To compile:

1. copy the cfacv/ directory wholesale under your home directory

2. cd ~/cfacv ; mkdir lib bin (if you haven't already)

3. cd src; make all

Note that this uses ``swig`` to build a shared-object library that the TcL interpreter in NAMD will load.  For this reason, it is a good idea if the TcL version that ``swig`` uses matches the TcL version that is "inside" NAMD.  If there is a mismatch, a good way to fix it is to compile NAMD from source using the system TcL header and libraries rather than the TcL the NAMD compilation instructions suggest.

In any run directory where you want to use `cfacv`, make sure the
   NAMD configuration file references the right path for `cfacv_tclforces.tcl`
   as being `$env(HOME)/cfacv/tcl/cfacv_tclforces.tcl` or wherever you chose to put 
   this repository.

### CAN I GET THIS TO WORK ON STAMPEDE2?

Yes!  Jim Phillips' most recent build of NAMD 2.12 (as of this writing, `/work/00288/tg455591/NAMD_LATEST_Linux-KNL-MPI-smp-Stampede`) cannot use CFACV because of a TcL incompatibility.  However, there is a version of `namd2` in my Stampede2 directory that _does_ work with CFACV, and it can be found at `/home1/00634/tg457991/namd/NAMD_2.12_Source/Linux-KNL-icc/namd2`.

### MORE DETAILS

CFACV interfaces with NAMD via a tclforces script called
`cfacv_tclforces.tcl` that initializes a workspace and then defines
the all-important `calcforces` procedure.  Procedures called from this
script are defined in cfacv.tcl and much of the implementation is
housed in cfacv.c.  Swig is used to interface cfacv.tcl and cfacv.c
to create a shared object library 'cfacv.so' loaded by
cfacv_tclforces.  (Note: you need `swig` installed to build the
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
beginning with cfacv_tclforces.tcl.  Issue `make cfacv.so` to build the
shared object library.

If you are only going to accelerate Cartesian variables, this is good
enough.  If you are going to accelerate dihedral angles, you have to
indicate that the restraint function in periodic with the additional
TAMD parameter {rf PERIODIC}.

Please contact me with any questions or problems.  This is not a professional
software release, so I fully anticipate there are bugs or other confusing
issues.

## ACKNOWLEDGMENTS

1. VMD and NAMD are products of the [Theoretical and Computational Biophysics Group at the NIH Center for Macromolecular Modeling and Bioinformatics at the University of Illinois at Urbana-Champaign](http://www.ks.uiuc.edu)

2. All codes and data in this repository have been made possible with partial support from NIH through grants AI084117, AI093248, GM115249, GM056550, and GM100472, the National Science Foundation through grants DMR-1207389 and MCB-1330205, and the US Army through grants W911NF-12-2-0022, W911-NF-13-1-0046, and W911NF-12-R-0011.

2017-2018, Cameron F Abrams

cfa22@drexel.edu

