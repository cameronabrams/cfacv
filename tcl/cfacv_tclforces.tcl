# tclforces file for CFACV
# Cameron F Abrams
# 2009-16

# check for base directory name variable;
# if not set, use default
if {![info exists CFACV_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(CFACV_BASEDIR)]} {
      set CFACV_BASEDIR $env(CFACV_BASEDIR)
  } else {
      set CFACV_BASEDIR ${HOME}/cfacv
  }
}

# will die if CFACV_BASEDIR is invalid
source ${CFACV_BASEDIR}/tcl/cfacv.tcl

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

# set the pseudorandom number generator seed
if {![info exists seed]} {
    set seed [exec /bin/date +%N | sed s/^0/1/]
    puts "CFACV) setting seed to $seed"
}

# declare and allocate data space
set ds [Tcl_NewDataSpace $nCntr $cvList $rList $seed]

# transfer the atom-center information (for now, just mass)
Tcl_DataSpace_SetCenterMasses $ds $masses

set first 1

# set output frequencies
if {[info exists TAMDof]} {
    set reportFreq $TAMDof
} else {
    set reportFreq 1
}

if {[info exists TAMDbinof]} {
    set binReportFreq $TAMDbinof
} else {
    set binReportFreq 1
}

if {[info exists restartINP]} {
    set first [Tcl_Reinitialize $ds $restartINP]
}

if {![info exists TAMDoutputlevel]} {
    set TAMDoutputlevel 3; # default output Z and Theta
}

set TAMDoutputFileFP 0
if {[info exists TAMDoutputFile]} {
    set TAMDoutputFileFP [my_fopen $TAMDoutputFile "w"]
    puts "CFACV) Opened TAMD output file $TAMDoutputFile"
} else {
    set TAMDoutputFileFP [my_fopen stdout "w"]
    puts "CFACV) TAMD output to stdout"
}

set TAMDbinOutputFileFP 0
if {[info exists TAMDbinOutputFile]} {
    set TAMDbinOutputFileFP [my_binfopen $TAMDbinOutputFile "w" $TAMDoutputlevel $ds]
    puts "CFACV) Binary TAMD output to $TAMDbinOutputFile"
}

# define the "calcforces" function which is called 
# at each MD timestep
proc calcforces { } {

    global ds
    global groups
    global first
    global reportFreq
    global binReportFreq
    global TAMDoutputlevel
    global TAMDoutputFileFP
    global TAMDbinOutputFile
    global TAMDbinOutputFileFP

    # load COM coordinates of requested atom groups into associative array
    loadcoords p

    # perform the update that transmits forces
    Tcl_UpdateDataSpace $ds p $groups $first [getstep]
    if {$first==1} { set first 0 }

    # report if requested
    if {[expr {[getstep]%$reportFreq == 0}]} {
	DataSpace_ReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDoutputFileFP
    }

    # report if requested
    if {[info exists TAMDbinOutputFile] && [expr {[getstep]%$binReportFreq == 0}]} {
	DataSpace_BinaryReportRestraints $ds [getstep] $TAMDoutputlevel $TAMDbinOutputFileFP
    }
}
