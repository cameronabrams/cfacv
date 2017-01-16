# tclforces file for CFACV implementing string method using NAMD's built-in replica implementation
# Cameron F Abrams
# 2014-2016

# below should be called by the main string method routine
# check for base directory name environment variable;
# if not set, use default
#if {![info exists CFACV_BASEDIR]} {
#  if {[info exists env(CFACV_BASEDIR)]} {
#    set CFACV_BASEDIR $env(CFACV_BASEDIR)
#  } else {
#    set HOME $env(HOME)
#    set CFACV_BASEDIR ${HOME}/cfacv
#  }
#}
#
# will die if CFACV_BASEDIR is invalid
#source ${CFACV_BASEDIR}/tcl/cfacv.tcl

cfacv_banner NAMD

# Check for necessary parameters set in main config. file
set tripped 0
foreach key {labelPDB cvINP} {
    if {![info exists $key]} {
	puts "CFACV) ERROR: you must set $key in the NAMD config. file."
	set tripped 1
    }
}
if {$tripped} {
    exit
}

set serArray {}; # must have for addgroup
set masses {}

# read the template PDB file that identifies subdomain memberships
set nCntr [read_centersPDB $labelPDB serArray masses atnames numatoms]
puts "CFACV) nCenters $nCntr  masses $masses"

# Set up the subdomains as "groups" for tclforces
set groups {}
for {set i 0} { $i < $nCntr } { incr i } {
    if {[info exists TAMDverbose]} {
	puts "addgroup $i :"
	puts "   [lindex $serArray $i]"
    }
    set gid [addgroup [lindex $serArray $i]]
    lappend groups $gid
}

# Set up list of CV's
set cvList {}
set nCV [read_cvs $cvINP cvList $nCntr]
puts "CFACV) nCV $nCV"
if {[info exists TAMDverbose]} {
    puts "CFACV) cvList: $cvList"
}

if {!$nCV} {
    puts "CFACV) ERROR: Perhaps you need to use mk_tPDB.tcl to generate the cv.inp file?"
}

# Set up list of restraints
set rList {}
# if the file listing explicit CV restraint matrix exists, use it
if {[info exists restrINP]} {
    set nR [read_restraints $restrINP $nCV rList]
    puts "CFACV) $restrINP : nRestraints $nR"
} else {
    set nR [create_single_cv_restraints $nCV rList $restrPARAMS]
    puts "CFACV) $nR single-cv restraints created:"
}

foreach r $rList {
    puts "CFACV) $r"
}

# set the pseudorandom number generator seed for this replica
# if it is not already set
if {![info exists seed]} {
    set seed [exec /bin/date +%N | sed s/^0/1/]
    puts "CFACV) setting seed to $seed"
}

# declare and allocate dataspace for this replica
set ds [Tcl_NewDataSpace $nCntr $cvList $rList $seed]

# transfer the atom-center information (for now, just mass)
Tcl_DataSpace_SetCenterMasses $ds $masses

# set up the metric tensor M
DataSpace_metricTensor_Setup $ds
DataSpace_metricTensor_Reset $ds

#set z {}
#set new_z 1 ; # signals that we need to set the Z's to the CV values

# define the "calcforces" function which is called 
# at each MD timestep
proc calcforces { } {
    global ds
    global groups
#    global reportFreq
#    global binReportFreq
#    global TAMDoutputlevel
#    global TAMDoutputFileFP
#    global TAMDbinOutputFile
#    global TAMDbinOutputFileFP
#    global SMPARAMS ; # associative array of string method parameters

    global new_z          ; # flag set to TRUE if this step is the first in an new iter (force and MT accum resets needed)

    global z        ; # z values to be communicated from master

#    global g        ; # dF/dz estimates to be sent back to master
#    global M        ; # metric tensor 
#    global steps_per_run  ; # number of steps per iteration; if [getstep]%nstepsperiter == 0, then this is the last step in the iter and a raw CV update is peformed
                          ; # to be fed back to the replica script

#    global k ;# spring constants
#    global th ;# theta's

    set timestep [getstep]

    # load COM coordinates of requested atom groups into associative array
    loadcoords p

    # given the z-values set by the master, determine the forces acting on the atoms,
    # and tally averages of dF/dz and M
    #puts stdout "New z? $new_z"
    #flush stdout
    Tcl_DataSpace_UpdateSM $ds p $groups $z $timestep $new_z
    # turn off the "new z" flag if it was on; it will be
    # turned back on by the main stringmethod tcl 
    # script when it updates the z's
    if { $new_z } {
       set new_z 0
    }
}
