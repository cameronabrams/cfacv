# Check the string-method machinery without NAMD: the image update
# z <- z - stepsize * M . dF/dz, reparameterization to equal arc length, and
# dual-image mode.
#
# Reparameterization is checked against an independent Tcl implementation of
# the same algorithm (piecewise-linear interpolation of the images to equal
# arc length, repeated until the spread of the inter-image distances falls
# below the tolerance), on a string that is curved and unevenly spaced, so
# that interpolation errors cannot cancel.
#
# Usage:  tclsh tests/string_method.tcl [CFACV_BASEDIR]
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

proc dist {a b} {
    set s 0.0
    foreach x $a y $b { set s [expr {$s + ($x-$y)**2}] }
    expr {sqrt($s)}
}

# A string of n images in a 2-D CV space.
proc new_sm {images {evolve_ends 1} {dual 0} {tol 1e-8} {maxiter 1000}} {
    set ni [llength $images]
    set sm [New_stringMethod_Dataspace $ni 2 0 0.0 $evolve_ends $dual]
    SMDataSpace_set_reparam_tol $sm $tol $maxiter
    set_images $sm $images
    return $sm
}
proc set_images {sm images} {
    set i 0
    foreach img $images {
        set z [SMDataSpace_image_z $sm $i]
        array_setitem $z 0 [lindex $img 0]
        array_setitem $z 1 [lindex $img 1]
        incr i
    }
}
proc get_images {sm ni} {
    set out {}
    for {set i 0} {$i < $ni} {incr i} {
        set z [SMDataSpace_image_z $sm $i]
        lappend out [list [array_getitem $z 0] [array_getitem $z 1]]
    }
    return $out
}
# dF/dz and M for every image
proc set_forces {sm ni g M} {
    for {set i 0} {$i < $ni} {incr i} {
        set ga [SMDataSpace_image_g $sm $i]
        for {set k 0} {$k < 2} {incr k} { array_setitem $ga $k [lindex $g $k] }
        set Ma [SMDataSpace_image_M $sm $i]
        for {set k 0} {$k < 4} {incr k} { array_setitem $Ma $k [lindex $M $k] }
    }
}

# Reference: one pass of linear interpolation to equal arc length.
proc ref_reparam_pass {imgs} {
    set n [llength $imgs]
    set L {0.0}
    for {set i 1} {$i < $n} {incr i} {
        lappend L [expr {[lindex $L end] + [dist [lindex $imgs $i] [lindex $imgs [expr {$i-1}]]]}]
    }
    set dL [expr {[lindex $L end] / ($n - 1)}]
    set out [list [lindex $imgs 0]]
    for {set i 1} {$i < $n-1} {incr i} {
        set s [expr {$i * $dL}]
        for {set k 1} {$k < $n} {incr k} {
            if {[lindex $L [expr {$k-1}]] < $s && $s <= [lindex $L $k]} break
        }
        set a [lindex $imgs [expr {$k-1}]]; set b [lindex $imgs $k]
        set f [expr {($s - [lindex $L [expr {$k-1}]]) / ([lindex $L $k] - [lindex $L [expr {$k-1}]])}]
        lappend out [lmap x $a y $b {expr {$x + $f*($y - $x)}}]
    }
    lappend out [lindex $imgs end]
    return $out
}
proc spacing_spread {imgs} {
    set n [llength $imgs]
    set sum 0.0; set sum2 0.0
    for {set i 1} {$i < $n} {incr i} {
        set d [dist [lindex $imgs [expr {$i-1}]] [lindex $imgs $i]]
        set sum [expr {$sum + $d}]; set sum2 [expr {$sum2 + $d*$d}]
    }
    expr {sqrt(max(0.0, ($sum2 - $sum*$sum/($n-1)) / ($n-1)))}
}
proc ref_reparam {imgs tol maxiter} {
    for {set it 0} {$it < $maxiter} {incr it} {
        set imgs [ref_reparam_pass $imgs]
        if {[spacing_spread $imgs] <= $tol} break
    }
    return $imgs
}

proc check_images {name got want tol} {
    set i 0
    foreach g $got w $want {
        check "$name image $i x" [lindex $g 0] [lindex $w 0] $tol
        check "$name image $i y" [lindex $g 1] [lindex $w 1] $tol
        incr i
    }
}

set PI [expr {acos(-1.0)}]
set I2 {1.0 0.0 0.0 1.0}
set zero {0.0 0.0}

# An unevenly spaced quarter circle: the angles bunch toward the start.
set n 9
set curve {}
for {set i 0} {$i < $n} {incr i} {
    set t [expr {($PI/2) * ($i/double($n-1))**2}]
    lappend curve [list [expr {cos($t)}] [expr {sin($t)}]]
}

# 1. Reparameterization alone (zero force), against the reference.
set sm [new_sm $curve]
set_forces $sm $n $zero $I2
SMDataSpace_MoveString $sm 0.0
set got [get_images $sm $n]
check_images "reparam" $got [ref_reparam $curve 1e-8 1000] 1e-9
check "reparam spacing spread" [spacing_spread $got] 0.0 1e-8

# 2. A second reparameterization, of the string the first one produced, must
#    also match the reference: the first call must not leave state behind
#    that changes the second.
SMDataSpace_MoveString $sm 0.0
check_images "reparam twice" [get_images $sm $n] [ref_reparam $got 1e-8 1000] 1e-12

# 3. The update: an evenly spaced straight string under a uniform force moves
#    rigidly by -stepsize * M . g (reparameterization then has nothing to do).
set line {}
for {set i 0} {$i < 5} {incr i} { lappend line [list [expr {0.25*$i}] 0.5] }
set M {2.0 0.5 0.5 1.0}
set g {0.3 -0.8}
set h 0.01
set dz [list [expr {-$h*(2.0*0.3 + 0.5*-0.8)}] [expr {-$h*(0.5*0.3 + 1.0*-0.8)}]]
set sm [new_sm $line]
set_forces $sm 5 $g $M
SMDataSpace_MoveString $sm $h
set want [lmap p $line {lmap x $p d $dz {expr {$x + $d}}}]
check_images "update" [get_images $sm 5] $want 1e-12

# 4. With fixed ends, the ends stay put and the interior moves.
set sm [new_sm $line 0]
set_forces $sm 5 $g $M
SMDataSpace_MoveString $sm $h
set got [get_images $sm 5]
check_images "fixed ends: ends" [list [lindex $got 0] [lindex $got end]] [list [lindex $line 0] [lindex $line end]] 1e-12
check "fixed ends: interior moved" [expr {[dist [lindex $got 2] [lindex $line 2]] > 1e-4}] 1 0

# 5. Dual-image mode: the first half of the images is the string, and each
#    image in the second half follows its partner in the first half.
set sm [new_sm [concat $curve $curve] 1 1]
set_forces $sm [expr {2*$n}] $zero $I2
SMDataSpace_MoveString $sm 0.0
set got [get_images $sm [expr {2*$n}]]
check_images "dual: reference half" [lrange $got 0 [expr {$n-1}]] [ref_reparam $curve 1e-8 1000] 1e-9
check_images "dual: partners" [lrange $got $n end] [lrange $got 0 [expr {$n-1}]] 1e-12

# 6. A string that converges on the last allowed reparameterization
#    iteration has converged, and must not abort the run. (Run in a child
#    interpreter, because the failure mode is exit().)
set child [format {
    load %s cfa_cvlibc
    set sm [New_stringMethod_Dataspace 5 2 0 0.0 1 0]
    SMDataSpace_set_reparam_tol $sm 1e-3 1
    for {set i 0} {$i < 5} {incr i} {
        set z [SMDataSpace_image_z $sm $i]
        array_setitem $z 0 [expr {0.25*$i}]
        array_setitem $z 1 0.0
    }
    SMDataSpace_MoveString $sm 0.0
    puts RETURNED
} [list [file join $CFACV_BASEDIR lib cfacv.so]]]
catch {exec [info nameofexecutable] << $child 2>@1} out
check "converged at maxiter returns" [string match *RETURNED* $out] 1 0

puts "string_method: $npass passed, $nfail failed"
exit [expr {$nfail > 0}]
