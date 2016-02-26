#
# Example VMD script for use in generating
# a CG label PDB; copy and edit as needed
#
# 2009-16, Cameron F Abrams
#
# Invoke as
#
# vmd -dispdev text -e mk_tPDB.tcl -args -pdb <base_pdb_file_name> -o <output_pdb_file_name> -cv <cvoption>
#
# The output defaults to "label.pdb".
#
# cvoption can be either "all-to-all" or "shortest-N", where N is an integer.  In the 
# case of "all-to-all", mk_tPDB.tcl will generate a cv input file "cv.inp" that
# declares all possible intercenter bonds.  In the caes of "shortest-N", it will
# generate cv.inp in which only the shortest N bonds are declared.
#
# You should change the value of CFACV_BASEDIR below to correspond to the
# location of cfacv/ in your account.
#

# check for base directory name environment variable;
# if not set, use default
if {[info exists env(CFACV_BASEDIR)]} {
    set CFACV_BASEDIR $env(CFACV_BASEDIR)
} else {
#   set CFACV_BASEDIR $env(HOME)/cfacv
    set CFACV_BASEDIR .. ; # for purposes of this example
}

source ${CFACV_BASEDIR}/tcl/cfacv.tcl
source ${CFACV_BASEDIR}/tcl/block.tcl


# some command-line argument parsing
set pdb        [getArg $argv "pdb" "p" {}]
set targetName [getArg $argv "targ" "o" "label.pdb"]
set cvOpt      [split [getArg $argv "cv" "cv" {}] ,]
set minResiduesPerCenter [getArg $argv "mrpc" "n" 20]

cfacv_banner $argv

if {![string length $pdb]} {
    puts "CFACV) ERROR: You must provide a base PDB using -p"
    exit
}

if {[file exists $pdb]} {
    set at_mid [initialize_base_template $pdb]
} else { 
    puts "CFACV) ERROR: base PDB $pdb does not exist."
    exit
}


set ids [NewIGenDataSpace 100 500]
set gds [NewGenDataSpace  100 500]

set p {}

##############################
## begin user-editable code

block_by_rgyr p $at_mid "protein" 5

# If you want to add a single mapping point for a selection, 
# use lappend, as below:
# lappend p [atomselect $at_mid "my_domain_string"]
#
## end user-editable code
##############################


# 3. Assign beta values
assign_center_betas $p 0

# 4. Generate the PDB
set all [atomselect $at_mid "all"]
$all writepdb $targetName
puts "CFACV) Generated $targetName."

# 5. Write information about each center to an output file called "centers.dat"
report_centers $p "centers.dat"

# 6.  If specified, generate a cv.inp file
if {[llength $cvOpt]} {
    set cvList [generate_cv $cvOpt $p]
    output_cv $cvList "cv.inp"
    puts "CFACV) Generated cv.inp"
}

puts "CFACV) mk_tPDB.tcl ends."

exit
