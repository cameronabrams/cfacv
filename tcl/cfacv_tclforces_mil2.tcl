# tclforces file for CFACV
# Cameron F Abrams
# 2009-13

# check for base directory name environment variable;
# if not set, use default
if {[info exists env(CFACV_BASEDIR)]} {
    set CFACV_BASEDIR $env(CFACV_BASEDIR)
} else {
    set HOME $env(HOME)
    set CFACV_BASEDIR ${HOME}/cfacv
}

# will die if CFACV_BASEDIR is invalid
source ${CFACV_BASEDIR}/cfacv.tcl

cfacv_banner NAMD

# Check for necessary parameters set in main config. file
set tripped 0
foreach key {labelPDB cvINP voronoi_centers_file home_cell} {
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
puts "CFACV) nCenters $nCntr  masses $masses numatoms $numatoms"

puts "CFACV/MIL) Adding all atoms to TCL-forces request"
for {set i 1} {$i <= $numatoms} {incr i} {
    addatom $i
}
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

# read in the voronoi cell centers; 
# each center is an (nCV)-component vector (a location in CV-space)
set vcList {}
set nVC [read_vcs $voronoi_centers_file vcList $nCV]
puts "CFACV) nVC $nVC"

# set the pseudorandom number generator seed
if {![info exists seed]} {
    set seed [exec /bin/date +%N | sed s/^0/1/]
    puts "CFACV) setting seed to $seed"
}


# declare, allocate, and populate data space
set ds [Tcl_NewDataSpace_mil $nCntr $cvList $vcList $home_cell $seed]
set first 1

#puts "developmental friendly dataspace dump"
DataSpace_ReportAll $ds
flush stdout
#puts "developmental friendly exit"
#exit

# set output frequencies
if {[info exists cfacv_of]} {
    set reportFreq $cfacv_of
} else {
    set reportFreq 1
}

if {[info exists cfacv_binof]} {
    set binReportFreq $cfacv_binof
} else {
    set binReportFreq 1
}

#if {[info exists restartINP]} {
#    set first [Tcl_Reinitialize $ds $restartINP]
#}

if {![info exists cfacv_outputlevel]} {
    set cfacv_outputlevel 3
}

set cfacv_outputFileFP 0
if {[info exists cfacv_outputFile]} {
    set cfacv_outputFileFP [my_fopen $cfacv_outputFile "w"]
    puts "CFACV) Opened CFACV output file $cfacv_outputFile"
} else {
    set cfacv_outputFileFP [my_fopen stdout "w"]
    puts "CFACV) CFACV output to stdout"
}

set cfacv_binOutputFileFP 0
if {[info exists cfacv_binOutputFile]} {
    set cfacv_binOutputFileFP [my_binfopen $cfacv_binOutputFile "w" $cfacv_outputlevel $ds]
    puts "CFACV) Binary CFACV output to $cfacv_binOutputFile"
}

if {![info exists write_config_at_hit]} {
    set write_config_at_hit 0 
}

if {![info exists die_on_sphere]} {
    set die_on_sphere 0
}

set violationID 0
set last_cell $home_cell
set thisviolation_stepcount 0

# define the "calcforces" function which is called 
# at each MD timestep
proc calcforces { } {

    global numatoms
    global ds
    global groups
    global first
    global reportFreq
    global binReportFreq
    global cfacv_outputlevel
    global cfacv_outputFileFP
    global cfacv_binOutputFile
    global cfacv_binOutputFileFP
    global home_cell
    global last_cell
    global reverse_postwaitsteps
    global outputname
    global atnames
    global write_config_at_hit
    global die_on_sphere
    global violationID
    global thisviolation_stepcount
    
    setboundaryflag 0

    set step [getstep]

    # load COM coordinates of requested atom groups into associative array
    loadcoords p
    
    # perform the update that detects whether we should reverse velocities

    set curr_cell [Tcl_UpdateDataSpace_mil $ds p $groups $first $step $home_cell]
    #print "CFACV/MIL) $step : curr_cell $curr_cell"

    # if this is the first time-step, make sure the cv-pt is in the home cell
    if {$first==1 && $curr_cell != $home_cell} {
	print "CFACV/MIL) ERROR:  initial placement of cv-pt is not in home cell."
	print "CFACV/MIL)         home_cell $home_cell curr_cell $curr_cell"
	exit
    }
    # puts "CFACV/MIL) DEBUG $step curr_cell $curr_cell"

    # check to see if we need a velocity reversal
    if { $curr_cell != $home_cell} {
	print "CFACV/MIL) $step : home $home_cell last $last_cell curr $curr_cell"
	if { $last_cell == $home_cell } {
	    set thisviolation_stepcount 0
	    incr violationID
	    print "CFACV/MIL) $step : violation $violationID : attempted passage from $home_cell to $curr_cell"
	} else {
	    incr thisviolation_stepcount
	}
	print "CFACV/MIL) $step : violation $violationID : stepcount $thisviolation_stepcount : revwt $reverse_postwaitsteps"
	# only if enough steps have elapsed since last reversal...
	if { ![expr $thisviolation_stepcount % $reverse_postwaitsteps]} {

	    print "CFACV/MIL) $step : violation $violationID : stepcount-in-violation $thisviolation_stepcount : reversing all velocities."
	    setboundaryflag 1
	    DataSpace_ReportCVs $ds $step
	    
	    # write a config if requested
	    if {$write_config_at_hit} {
		set fp [open "${outputname}_${home_cell}_${curr_cell}_${step}.xyz" "w"]
		puts $fp "$numatoms"
		puts $fp " generated by cfacv/mil/tclforces at step $step"
		for {set i 1} {$i <= $numatoms} {incr i} {
		    puts $fp " [lindex $atnames [expr $i-1]] $p($i)"
		}
		close $fp
		print "CFACV/MIL) Wrote ${outputname}_${home_cell}_${curr_cell}_${step}.xyz"
	    }
	}
	if { $curr_cell == -1 && $die_on_sphere} {
	    print "CFACV/MIL)  CV is outside sphere"
	    print "CFACV/MIL)  Exit."
	    exit
	}
    }
    set last_cell $curr_cell
    if {$first==1} { set first 0 }

}
