# Check bin/FEPFromForces, which integrates the mean force along a string
# into a free-energy profile, against free energies known in closed form.
#
# The input force is f = -grad F at each image. Along a straight chord the
# force of a quadratic F is linear, so the trapezoid rule on the chords is
# exact for a quadratic F on any string, curved or not; so is the arc length,
# the sum of the chord lengths. A second string crosses the +/-pi seam of a
# periodic CV, which -dihed must wrap.
#
# Usage:  tclsh tests/fep_from_forces.tcl [CFACV_BASEDIR]
# Exits 0 when every check passes, 1 otherwise.

set CFACV_BASEDIR [expr {$argc > 0 ? [lindex $argv 0] : [file normalize [file join [file dirname [info script]] ..]]}]
set exe [file join $CFACV_BASEDIR bin FEPFromForces]
if {![file executable $exe]} {
    set exe [file join $CFACV_BASEDIR src FEPFromForces]
}

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

set tmp [file join [expr {[info exists ::env(TMPDIR)] ? $::env(TMPDIR) : "/tmp"}] fep_from_forces.[pid]]
file mkdir $tmp

# Run FEPFromForces on a list of {z f} image records; return the output
# rows as lists {j s F}.
proc run_fep {images args} {
    global exe tmp
    set in [file join $tmp in.dat]
    set out [file join $tmp out.dat]
    set fp [open $in w]
    puts $fp "# z f"
    foreach img $images { puts $fp [join [concat {*}$img] " "] }
    close $fp
    exec $exe -i $in -o $out -dim [llength [lindex $images 0 0]] {*}$args
    set fp [open $out r]
    set rows {}
    foreach line [split [string trim [read $fp]] \n] { lappend rows [lrange $line 0 2] }
    close $fp
    return $rows
}

# 1. F(x,y) = 3x^2 + 2xy - y^2 on a curved, unevenly spaced string.
proc F {x y} { expr {3*$x*$x + 2*$x*$y - $y*$y} }
set images {}
set pts {}
for {set j 0} {$j < 11} {incr j} {
    set t [expr {($j/10.0)**1.5}]
    set x [expr {cos(2.0*$t)}]
    set y [expr {sin(2.0*$t) + 0.3*$t}]
    lappend pts [list $x $y]
    # f = -grad F
    lappend images [list [list $x $y] [list [expr {-(6*$x + 2*$y)}] [expr {-(2*$x - 2*$y)}]]]
}
set rows [run_fep $images]
check "quadratic: number of images" [llength $rows] 11 0
set s 0.0
for {set j 0} {$j < 11} {incr j} {
    lassign [lindex $pts $j] x y
    lassign [lindex $pts 0] x0 y0
    if {$j > 0} {
        lassign [lindex $pts [expr {$j-1}]] xp yp
        set s [expr {$s + hypot($x-$xp, $y-$yp)}]
    }
    lassign [lindex $rows $j] jj sj Fj
    check "quadratic: image index $j" $jj $j 0
    check "quadratic: F at image $j" $Fj [expr {[F $x $y] - [F $x0 $y0]}] 1e-4
    check "quadratic: arc length at image $j" $sj $s 1e-4
}

# 2. A periodic CV across the seam: phi from 2.9 rad up through pi to
#    -2.9 rad, with F = c * phi (unwrapped), so f = -c.
set PI [expr {acos(-1.0)}]
set c 2.5
set images {}
for {set j 0} {$j < 9} {incr j} {
    set phi [expr {2.9 + 0.06*$j}]
    set wrapped [expr {$phi > $PI ? $phi - 2*$PI : $phi}]
    lappend images [list [list $wrapped] [list [expr {-$c}]]]
}
set rows [run_fep $images -dihed]
for {set j 0} {$j < 9} {incr j} {
    lassign [lindex $rows $j] jj sj Fj
    check "seam: F at image $j" $Fj [expr {$c*0.06*$j}] 1e-4
    check "seam: arc length at image $j" $sj [expr {0.06*$j}] 1e-4
}

# 3. Fewer than two images is an error, not a profile.
set images [list [list {0.0 0.0} {1.0 1.0}]]
check "one image is rejected" [catch {run_fep $images}] 1 0

file delete -force $tmp
puts "fep_from_forces: $npass passed, $nfail failed"
exit [expr {$nfail > 0}]
