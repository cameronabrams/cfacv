# Check targeted MD (TMD): the restraint target must stay put before step
# t0, move by equal increments, reach the TMD target exactly at step t1, and
# stay put after it -- for harmonic and for periodic restraints alike.
#
# Usage:  tclsh tests/tmd.tcl [CFACV_BASEDIR]
# Exits 0 when every check passes, 1 otherwise.

set CFACV_BASEDIR [expr {$argc > 0 ? [lindex $argv 0] : [file normalize [file join [file dirname [info script]] ..]]}]
load [file join $CFACV_BASEDIR lib cfacv.so] cfa_cvlibc

set nfail 0
set npass 0
proc check {label got want tol} {
    global nfail npass
    if {abs($got - $want) > $tol} {
        puts [format "FAIL %-50s got % .10f want % .10f" $label $got $want]
        incr nfail
    } else {
        incr npass
    }
}
proc intarray {l} {
    set a [new_arrayint [llength $l]]
    set i 0
    foreach x $l { arrayint_setitem $a $i $x; incr i }
    return $a
}
proc dblarray {l} {
    set a [new_array [llength $l]]
    set i 0
    foreach x $l { array_setitem $a $i $x; incr i }
    return $a
}

# One CV over n point-mass centers, one TMD restraint on it.
proc make_tmd {cvtype pos rf zmin zmax target t0 t1} {
    set n [llength $pos]
    set ds [NewDataSpace $n 1 1 12345]
    set ind {}
    for {set i 0} {$i < $n} {incr i} { DataSpace_AddCenter_MassOnly $ds 1.0; lappend ind $i }
    DataSpace_AddCV $ds $cvtype $n [intarray $ind]
    DataSpace_AddRestr $ds 10.0 0.0 1 [dblarray {1.0}] $rf $zmin $zmax
    DataSpace_AddTmdOpt $ds 0 $target $t0 $t1
    return $ds
}

# Run steps 0..nsteps-1 with the centers held fixed; return the target z in
# effect at each step (the z that step's force was computed against).
proc run_tmd {ds pos nsteps} {
    set zs {}
    for {set step 0} {$step < $nsteps} {incr step} {
        set i 0
        foreach p $pos {
            set r [DataSpace_centerPos $ds $i]
            for {set d 0} {$d < 3} {incr d} { array_setitem $r $d [lindex $p $d] }
            incr i
        }
        DataSpace_ComputeCVs $ds
        # z before this step's update is the target this step's force used
        if {$step > 0} { lappend zs [array_getitem [DataSpace_z $ds] 0] }
        DataSpace_RestrainingForces $ds [expr {$step == 0}] $step
    }
    lappend zs [array_getitem [DataSpace_z $ds] 0]
    return $zs
}

# For a periodic restraint, distances are wrapped into [-pi, pi]: the target
# goes the short way round, and values are compared modulo 2 pi.
proc wrap {x} {
    set pi [expr {acos(-1.0)}]
    while {$x > $pi} { set x [expr {$x - 2*$pi}] }
    while {$x < -$pi} { set x [expr {$x + 2*$pi}] }
    return $x
}
proc check_schedule {name zs z0 target t0 t1 {periodic 0}} {
    # zs[s] is the target in effect at step s+1 (after s+1 updates)
    set dist [expr {$target - $z0}]
    if {$periodic} { set dist [wrap $dist] }
    set inc [expr {$dist / double($t1 - $t0)}]
    for {set s 0} {$s < [llength $zs]} {incr s} {
        set step [expr {$s + 1}]
        if {$step <= $t0} {
            set want $z0
        } elseif {$step >= $t1} {
            set want $target
        } else {
            set want [expr {$z0 + ($step - $t0) * $inc}]
        }
        set got [lindex $zs $s]
        if {$periodic} {
            check "$name target in domain at step $step" [expr {abs($got) <= acos(-1.0)}] 1 0
            set got [expr {$want + [wrap [expr {$got - $want}]]}]
        }
        check "$name target at step $step" $got $want 1e-9
    }
}

set PI [expr {acos(-1.0)}]
set t0 5
set t1 15
set nsteps 25

# harmonic: a bond from 1.5 A to 2.5 A
set pos {{0.0 0.0 0.0} {1.5 0.0 0.0}}
set ds [make_tmd BOND $pos HARMONIC 0 0 2.5 $t0 $t1]
check_schedule BOND [run_tmd $ds $pos $nsteps] 1.5 2.5 $t0 $t1

# periodic: a dihedral of +60 degrees driven to +120 degrees
set c [expr {cos($PI/3)}]; set s [expr {sin($PI/3)}]
set pos [list {1.0 0.0 1.0} {0.0 0.0 1.0} {0.0 0.0 0.0} [list $c $s 0.0]]
# read the starting dihedral straight from the CV, not through a restraint
set ds [NewDataSpace 4 1 0 12345]
for {set i 0} {$i < 4} {incr i} { DataSpace_AddCenter_MassOnly $ds 1.0 }
DataSpace_AddCV $ds DIHED 4 [intarray {0 1 2 3}]
set i 0
foreach p $pos {
    set r [DataSpace_centerPos $ds $i]
    for {set d 0} {$d < 3} {incr d} { array_setitem $r $d [lindex $p $d] }
    incr i
}
DataSpace_ComputeCVs $ds
set res [new_array 1]
DataSpace_WriteCVtoArray $ds [intarray {1}] $res
set z0 [array_getitem $res 0]
check "DIHED starting value is +/-60 degrees" [expr {abs($z0)}] [expr {$PI/3}] 1e-9
set target [expr {2*$z0}]
set ds [make_tmd DIHED $pos PERIODIC -$PI $PI $target $t0 $t1]
check_schedule DIHED [run_tmd $ds $pos $nsteps] $z0 $target $t0 $t1 1

# periodic across the seam: +150 degrees to -150 degrees is 60 degrees the
# short way, through +/-180, not 300 degrees the long way
set a [expr {-5*$PI/6}]
set pos [list {1.0 0.0 1.0} {0.0 0.0 1.0} {0.0 0.0 0.0} [list [expr {cos($a)}] [expr {sin($a)}] 0.0]]
set ds [NewDataSpace 4 1 0 12345]
for {set i 0} {$i < 4} {incr i} { DataSpace_AddCenter_MassOnly $ds 1.0 }
DataSpace_AddCV $ds DIHED 4 [intarray {0 1 2 3}]
set i 0
foreach p $pos {
    set r [DataSpace_centerPos $ds $i]
    for {set d 0} {$d < 3} {incr d} { array_setitem $r $d [lindex $p $d] }
    incr i
}
DataSpace_ComputeCVs $ds
set res [new_array 1]
DataSpace_WriteCVtoArray $ds [intarray {1}] $res
set z0 [array_getitem $res 0]
check "DIHED seam starting value is +150 degrees" $z0 [expr {5*$PI/6}] 1e-9
set target [expr {-5*$PI/6}]
set ds [make_tmd DIHED $pos PERIODIC -$PI $PI $target $t0 $t1]
check_schedule "DIHED seam" [run_tmd $ds $pos $nsteps] $z0 $target $t0 $t1 1

puts "tmd: $npass passed, $nfail failed"
exit [expr {$nfail > 0}]
