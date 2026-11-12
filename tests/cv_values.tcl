# Check CFACV's collective variables without NAMD.
#
# For each CV type this loads lib/cfacv.so into a plain tclsh, places a few
# point-mass centers at known positions, and checks two things:
#
#   1. the CV value against one computed independently here, in Tcl;
#   2. the restraining force on every center, as returned through the
#      centers' position arrays exactly as tclforces would receive it,
#      against -dU/dx from a central finite difference of the restraint
#      energy.
#
# The second check is the one that matters for dynamics: it catches a wrong
# gradient, a sign error or a wrong chain-rule factor even when the value is
# right.
#
# Usage:  tclsh tests/cv_values.tcl [CFACV_BASEDIR]
# Exits 0 when every check passes, 1 otherwise.

set CFACV_BASEDIR [expr {$argc > 0 ? [lindex $argv 0] : [file normalize [file join [file dirname [info script]] ..]]}]
load [file join $CFACV_BASEDIR lib cfacv.so] cfa_cvlibc

set PI [expr {acos(-1.0)}]
set nfail 0
set npass 0

proc check {label got want tol} {
    global nfail npass
    set err [expr {abs($got - $want)}]
    set scale [expr {max(1.0, abs($want))}]
    if {$err > $tol * $scale} {
        puts [format "FAIL %-40s got % .10f want % .10f" $label $got $want]
        incr nfail
    } else {
        incr npass
    }
}

proc dblarray {l} {
    set a [new_array [llength $l]]
    set i 0
    foreach x $l { array_setitem $a $i $x; incr i }
    return $a
}

proc intarray {l} {
    set a [new_arrayint [llength $l]]
    set i 0
    foreach x $l { arrayint_setitem $a $i $x; incr i }
    return $a
}

# A DataSpace with the given center masses, one CV over all of them, and one
# restraint on that CV.
proc make_ds {masses cvtype k z rftype zmin zmax} {
    set n [llength $masses]
    set ds [NewDataSpace $n 1 1 12345]
    foreach m $masses { DataSpace_AddCenter_MassOnly $ds $m }
    set ind {}
    for {set i 0} {$i < $n} {incr i} { lappend ind $i }
    DataSpace_AddCV $ds $cvtype $n [intarray $ind]
    DataSpace_AddRestr $ds $k $z 1 [dblarray {1.0}] $rftype $zmin $zmax
    return $ds
}

proc set_positions {ds pos} {
    set i 0
    foreach p $pos {
        set r [DataSpace_centerPos $ds $i]
        for {set d 0} {$d < 3} {incr d} { array_setitem $r $d [lindex $p $d] }
        incr i
    }
}

# Evaluate at pos; returns {cv-value energy forces}.  RestrainingForces
# overwrites the position arrays with forces, so positions are set each time.
proc evaluate {ds pos} {
    set_positions $ds $pos
    DataSpace_ComputeCVs $ds
    set active [intarray {1}]
    set res [new_array 1]
    DataSpace_WriteCVtoArray $ds $active $res
    set v [array_getitem $res 0]
    DataSpace_RestrainingForces $ds 0 0
    set u [DataSpace_RestraintEnergy $ds]
    set f {}
    for {set i 0} {$i < [llength $pos]} {incr i} {
        set r [DataSpace_centerPos $ds $i]
        lappend f [list [array_getitem $r 0] [array_getitem $r 1] [array_getitem $r 2]]
    }
    return [list $v $u $f]
}

proc check_cv {name ds pos want_value} {
    lassign [evaluate $ds $pos] v u f
    check "$name value" $v $want_value 1e-10
    set h 1e-6
    for {set i 0} {$i < [llength $pos]} {incr i} {
        for {set d 0} {$d < 3} {incr d} {
            set pp $pos; set pm $pos
            lset pp $i $d [expr {[lindex $pos $i $d] + $h}]
            lset pm $i $d [expr {[lindex $pos $i $d] - $h}]
            set up [lindex [evaluate $ds $pp] 1]
            set um [lindex [evaluate $ds $pm] 1]
            set fd [expr {-($up - $um) / (2.0*$h)}]
            check "$name force center $i dim $d" [lindex $f $i $d] $fd 1e-5
        }
    }
}

# --- independent reference implementations ---------------------------------

proc vsub {a b} { lmap x $a y $b {expr {$x - $y}} }
proc vdot {a b} { set s 0.0; foreach x $a y $b { set s [expr {$s + $x*$y}] }; return $s }
proc vnorm {a} { expr {sqrt([vdot $a $a])} }
proc vcross {a b} {
    lassign $a a0 a1 a2; lassign $b b0 b1 b2
    list [expr {$a1*$b2 - $a2*$b1}] [expr {$a2*$b0 - $a0*$b2}] [expr {$a0*$b1 - $a1*$b0}]
}
proc ref_angle {p0 p1 p2} {
    set a [vsub $p0 $p1]; set b [vsub $p2 $p1]
    expr {acos([vdot $a $b] / ([vnorm $a] * [vnorm $b]))}
}
# IUPAC convention (trans = +/-pi), as VMD's "measure dihed" reports it;
# Blondel & Karplus, J. Comput. Chem. 17, 1132 (1996).
proc ref_dihed {p0 p1 p2 p3} {
    set b0 [vsub $p1 $p0]; set b1 [vsub $p2 $p1]; set b2 [vsub $p3 $p2]
    expr {atan2([vnorm $b1] * [vdot $b0 [vcross $b1 $b2]],
                [vdot [vcross $b0 $b1] [vcross $b1 $b2]])}
}
proc ref_rgyr {pos masses} {
    set mt 0.0; set c {0.0 0.0 0.0}
    foreach p $pos m $masses {
        set mt [expr {$mt + $m}]
        set c [lmap ci $c pi $p {expr {$ci + $m*$pi}}]
    }
    set c [lmap ci $c {expr {$ci / $mt}}]
    set s 0.0
    foreach p $pos m $masses { set d [vsub $p $c]; set s [expr {$s + $m*[vdot $d $d]}] }
    expr {sqrt($s / $mt)}
}

# --- the checks ---------------------------------------------------------------

# Geometries are deliberately generic: no right angles, no coplanar points, no
# coordinate that happens to be zero, so a dropped term cannot cancel out.
set p0 {0.31 -1.17 0.42}
set p1 {1.52 -0.23 0.11}
set p2 {2.04 1.06 0.87}
set p3 {3.41 1.35 0.02}

set pos [list $p0 $p1]
check_cv BOND [make_ds {1.0 2.0} BOND 10.0 1.0 HARMONIC 0 0] $pos [vnorm [vsub $p0 $p1]]

set pos [list $p0 $p1 $p2]
check_cv ANGLE [make_ds {1.0 2.0 3.0} ANGLE 10.0 1.5 HARMONIC 0 0] $pos [ref_angle $p0 $p1 $p2]

set pos [list $p0 $p1 $p2 $p3]
set want [ref_dihed $p0 $p1 $p2 $p3]
check_cv DIHED [make_ds {1.0 2.0 3.0 4.0} DIHED 10.0 0.5 PERIODIC -$PI $PI] $pos $want
# A target across the +/-pi seam from the value: a periodic restraint must take
# the short way round, so the energy stays smooth and the forces still match.
set zseam [expr {$want > 0 ? $want - 2*$PI + 0.3 : $want + 2*$PI - 0.3}]
check_cv "DIHED across seam" [make_ds {1.0 2.0 3.0 4.0} DIHED 10.0 $zseam PERIODIC -$PI $PI] $pos $want

set masses {1.0 2.0 3.0 4.0}
check_cv RGYR [make_ds $masses RGYR 10.0 0.5 HARMONIC 0 0] $pos [ref_rgyr $pos $masses]

foreach typ {CARTESIAN_X CARTESIAN_Y CARTESIAN_Z} d {0 1 2} {
    check_cv $typ [make_ds {1.0} $typ 10.0 0.0 HARMONIC 0 0] [list $p0] [lindex $p0 $d]
}

puts "cv_values: $npass passed, $nfail failed"
exit [expr {$nfail > 0}]
