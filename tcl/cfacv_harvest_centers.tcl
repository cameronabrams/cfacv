# CFACV -- Collective Variables for NAMD
# Single-Sweep Free Energy reconstruction: Harvesting Centers
#
# This TcL script harvests centers from TAMD trajectories
# to build an irregular mesh in CV space
# 
# Invoke within VMD; e.g.,
#
# vmd -dispdev text -e cfacv_harvest_centers.tcl -args -psf XXX.psf -dcd 1.dcd,2.dcd -cv cv.inp -l label.pdb -restr restr.inp -stip DIST
#
# In the above, XXX.psf is the base PSF file, 1.dcd and 2.dcd are TAMD trajectories,
# cv.inp, label.pdb, and restr.inp are the CFACV input files used in
# the TAMD simulations (they define the CV's), and DIST is the minimum
# distance between harvested centers.
#
# Shared object libraries required:
#
# $CFACV_BASEDIR/lib/cfacv.so
# $CFACV_BASEDIR/lib/libgenericdataspace.so
#
# Cameron F Abrams
# 2010-2016

# check for base directory name variable;
# if not set, use default
if {![info exists CFACV_BASEDIR]} {
  # see if user set an environment variable
  if {[info exists env(CFACV_BASEDIR)]} {
      set CFACV_BASEDIR $env(CFACV_BASEDIR)
  } else {
      set CFACV_BASEDIR $env(HOME)/research/cfacv
  }
}

# will die if CFACV_BASEDIR is invalid
source ${CFACV_BASEDIR}/tcl/cfacv.tcl
source ${CFACV_BASEDIR}/tcl/genericdataspace.tcl

# some default control variable values
# write the full coordinates for each center (yeah, you'd better do that)
set writeFrames 1
# make a coarse-grained PDB file (not really necessary)
set makeCentersPDB 0

# by default we expect there to be a cv.inp; a restr.inp is not necessary
# if there is as usual only one restraint applied per CV
#set cvINP    cv.inp
#set restrINP restr.inp
# note, I think my use of getArg below sets the defaults already

# happy birthday to me
set seed 22772

# some argument handling; getArg defined in cfacv.tcl
set psf      [getArg $argv "psf"    "psf"   {}   ]
set pdb      [getArg $argv "pdb"    "pdb"   {} ]
set nPDB [llength $pdb]
set dcd      [getArg $argv "dcd"    "dcd"   {}   ]
set nDCD [llength $dcd]
set coor     [getArg $argv "coor"    "coor"   {}   ]
set nCOOR [llength $coor]
set labelPDB [getArg $argv "label"  "l"   label.pdb]
set cvINP [getArg $argv "cv"  "cv"   cv.inp]
set restrINP [getArg $argv "restr" "restr" restr.inp] 
set step     [getArg $argv "step" "step" 1]
# by default, put frames coor files in the current directory
set centerdir [getArg $argv "centerdir" "cd" .]
# stip is the minimum distance between any two centers, in angstroms
set stip     [getArg $argv "stip" "stip"  2.5]
# prostip is the maximum distance of any center to any protein atom, in angstroms
# since this only makes sense for cartesian CV's, let's not use it at all for now
#set prostip  [getArg $argv "prostip" "prostip" 5.0]

# overwrite any arg values with a file that can set variables 
if {[file exists cfacv_harvest_centers_local.tcl]} {
    source cfacv_harvest_centers_local.tcl
}

cfacv_banner $argv

# a procedure that can output the CG PDB file
proc make_PDB { cv typ name } {
    
    set nCenters [llength $cv]
    puts "CFACV) Writing [llength $cv] coords to \[$name\]."
#    puts "CFACV) DB: cv $cv"
#    puts "CFACV) DB: typ $typ"
    set fo [open $name "w"]
    puts $fo "REMARK CFACV 2010"
    for {set i 0} { $i < $nCenters } { incr i } {
	set thiscm [lindex $cv $i]
#	puts "CFACV) DB cm $i $thiscm"
	set x [format "% 8.3f" [lindex $thiscm 0]]
	set y [format "% 8.3f" [lindex $thiscm 1]]
	set z [format "% 8.3f" [lindex $thiscm 2]]
#	if {[expr $i < 1000]} {
#	    set at "C[format %03i $i]"
#	} else {
	    set at "C[format %03i [expr $i % 1000]]"
#	}
	set ind [format "% 7i" $i]
	set bet [format "%5.2f" [lindex $typ $i]]
	puts $fo "ATOM${ind} $at CGX A   1    ${x}${y}${z}  1.00 $bet      A"
    }
    puts $fo "END"
    close $fo
}


# check for minimal required input files
foreach filename {psf labelPDB cvINP } {
    if {[string length [subst $$filename]] > 0 && ![file exists [subst $$filename]]} {
	puts "CFACV) ERROR: required file [subst $$filename] not found"
	exit
    }
}

# Need EITHER a PSF or a PDB; create the new VMD molecule from that
if {[string length $psf] > 0} {
    mol new $psf
} elseif {[string length [lindex $pdb 0]] > 0} {
    mol new $pdb
} else {
    puts "CFACV) ERROR:  cfacv_harvest_centers.tcl needs minimally either a PSF or a PDB file."
    exit
}

# build the list of input coordinate files based on order stipulated at runtime
set coords [concat $pdb $coor $dcd]
puts "CFACV) INFO: [llength $coords] coord-files:"
set i 0
foreach c $coords {
    puts "CFACV)   $i \[$c\]"
    incr i
}

set at_mid [molinfo top get id]

set serArray {}
set masses {}

# read the label.pdb; i.e., the template PDB file that identifies 
# group memberships (these are referred to as "centers" unfortunately.
# Here, we mostly reserve "center" to refer to a CV instance that is
# harvested from the TAMD sweep.
set nGroups [read_centersPDB $labelPDB serArray masses atnames numatoms]
set groupSel [center_selections $at_mid $serArray]
set i 0
foreach sel $groupSel {
    puts "CFACV) INFO: group $i has [$sel num] atoms"
    incr i
}
puts "CFACV) INFO: nGroups $nGroups  masses $masses"
# read the cv.inp file
set cvList {}
set nCV [read_cvs $cvINP cvList $nGroups]
puts "CFACV) INFO: nCV $nCV cvList: $cvList"

# check to see if any cv is periodic and if so,
# set a flag for it and establish the domain
# Currently, only DIHED-type CV's are supported here.
set periodic {}
set domain {}
for {set i 0} {$i < $nCV} {incr i} {
   set cv [lindex $cvList $i]
   set cvtyp [lindex $cv 0]
   if {$cvtyp == "DIHED"} {
     lappend periodic 1
     lappend domain [list [expr -2*acos(0.0)] [expr 2*acos(0.0)]]
   } else {
     lappend periodic 0
     lappend domain [list 0 0]
   }
}
puts "Periodic: $periodic"
puts "Domain: $domain"
# Set up list of restraints; could either be stipulated in restr.inp
# or assumed to be one per CV with parameters in the list restrPARAMS

set rList {}
# if the file listing explicit CV restraint matrix exists, use it
if {[string length $restrINP] > 0} {
    set nR [read_restraints $restrINP $nCV rList]
    puts "CFACV) $restrINP : nRestraints $nR"
} else {
    set nR [create_single_cv_restraints $nCV rList $restrPARAMS]
    puts "CFACV) $nR single-cv restraints created:"
}

# only count CV's for which TAMD friction is NOT set to 0
# this is to allow us to sample restrained CV spaces
set activeCV {}
foreach r $rList {
    puts "CFACV) INFO: [restr_getopt $r TAMDgamma TAMDfriction 0] [cvc_getcvi [restr_getcvc $r]]"
    lappend activeCV [expr {([restr_getopt $r TAMDgamma TAMDfriction 0])!=0.0}]
}
puts "CFACV) INFO: activeCV: $activeCV"

# I don't know what this is for
#set activeSel {}
#set s 0
#for {set i 0} {$i < [llength $activeCV] } { incr i 3 } {
#    if {[lindex $activeCV $i] && [lindex $activeCV [expr $i + 1]] && [lindex $activeCV [expr $i + 2]]} {
#	lappend activeSel [lindex $cntrSel $s]
#    }
##    incr s
#}

#puts "CFACV) INFO: number of active group selections: [llength $activeSel]"

set ds [Tcl_NewDataSpace $nGroups $cvList $rList $seed]
set ids [NewIGenDataSpace $nCV $nR]
set rds [NewGenDataSpace $nCV $nR]

set all [atomselect top "all"]
set pro [atomselect top "protein"]

set cnf 0     ;# cumulative frame count
set ncf 0     ;# coord file index
set nf_arr {} ;# number of frames in each coordinate file

set nfp 0 ;# total number of frames processed
set nc  0 ;# number of centers accepted
set member_CV {} ;# CV's accepted
set member_typ {} ;# index of each member's coord file
set cf_membercount {}  ;# count of members extracted from each coord file

set all_CV {} ;# CV's from all frames considered
set all_typ {} ;# index of each frame's coord file

foreach cf $coords {

    puts "CFACV) Working on $ncf \[$cf\]..."

    mol addfile $cf waitfor all autobonds off
    set inf [molinfo top get numframes]
    set cnf [expr $cnf + $inf]
    lappend nf_arr $inf
    set this_acc_count 0

    for { set i 0 } { $i < $inf } { incr i } {
	$pro frame $i
	$all frame $i
	Tcl_ObserveDataSpace $ds $groupSel $i
	Tcl_IGenDataSpace_Write $ids 0 $activeCV
	DataSpace_WriteCVtoArray $ds [IGenDataSpace_getAddr $ids 0] [GenDataSpace_getAddr $rds 0]
	set current_CV [ArrayToList [GenDataSpace_getAddr $rds 0] $nCV]
	lappend all_CV $current_CV
	lappend all_typ $ncf

	#puts "CFACV) DEBUG: frame $i current_CV $current_CV"
	
	set tripped 0
	for { set j 0 } { $j < $nc && !$tripped } { incr j } {
            # use minimum image convention
	    set crit [veclength [vecsub_mic $current_CV [lindex $member_CV $j] $periodic $domain]]
	   # puts "CFACV) DEBUG: $cf $i -> $j/$nc = $crit"
	    if { $crit < $stip } {
		set tripped 1 
	    }
	}
	if { !$tripped || !$nc } {  ;# this center can be accepted
	    #if { [within_range $activeSel $pro $prostip] } {
		lappend member_CV $current_CV
		lappend fi $nfp
		lappend member_typ $ncf
		incr this_acc_count
		if { $writeFrames } {
		    $all writenamdbin "${centerdir}/center[format %05i $nc].coor"
		    puts "CFACV) Frame $i in $cf is accepted center $nc => ${centerdir}/center[format %05i $nc].coor"
		}
		incr nc
	    #}
	}
	incr nfp
    }
    lappend cf_membercount $this_acc_count
    ;# Delete these frames
    animate delete beg 0 end [expr $inf - 1] 0
    puts "CFACV) Done with \[$cf\]; [llength $cf_membercount] so far."
    incr ncf
}
puts "CFACV) INFO: processed $nfp frames and harvested $nc centers."

if { $makeCentersPDB } { 

    make_PDB $member_CV $member_typ "centers.pdb"
    make_PDB $all_CV $all_typ "all_cvs.pdb"

}

set fp [open "mesh.dat" "w"]
foreach v $member_CV {
   for {set i 0} {$i < [llength $v]} {incr i} {
     puts -nonewline $fp "[lindex $v $i] "
   }
   puts $fp "" 
}
close $fp

foreach d $coords t $cf_membercount {
    puts "\[$d\] contributed $t centers."
}



quit


