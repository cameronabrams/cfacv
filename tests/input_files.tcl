# Check that cfacv reads its three input files -- the label PDB, cv.inp and
# restr.inp -- exactly as docs/source/reference.md ("Input files") specifies.
# Those formats are a public interface: other programs write them. A change
# that makes this test fail is an incompatible change to that interface.
#
# Usage:  tclsh tests/input_files.tcl [CFACV_BASEDIR]
# Exits 0 when every check passes, 1 otherwise.

set CFACV_BASEDIR [expr {$argc > 0 ? [lindex $argv 0] : [file normalize [file join [file dirname [info script]] ..]]}]
# the Tcl library is in the checkout under test even when the build is not
source [file join [file normalize [file join [file dirname [info script]] ..]] tcl cfacv.tcl]

set nfail 0
set npass 0
proc check {label got want} {
    global nfail npass
    if {$got ne $want} {
        puts "FAIL $label: got {$got} want {$want}"
        incr nfail
    } else {
        incr npass
    }
}
proc checknum {label got want tol} {
    global nfail npass
    if {abs($got - $want) > $tol} {
        puts "FAIL $label: got $got want $want"
        incr nfail
    } else {
        incr npass
    }
}

set tmp [file join [expr {[info exists ::env(TMPDIR)] ? $::env(TMPDIR) : "/tmp"}] cfacv_input_files.[pid]]
file mkdir $tmp
proc writefile {name text} {
    global tmp
    set f [file join $tmp $name]
    set fp [open $f w]
    puts -nonewline $fp $text
    close $fp
    return $f
}

# One fixed-column PDB record. occ and beta are passed as the literal
# 6-character strings to put in columns 55-60 and 61-66.
proc pdbline {rec serial name resname resid occ beta} {
    # columns: record 1-6, serial 7-11, name 13-16, resname 18-21, chain 22,
    # resid 23-26, x y z 31-54, occupancy 55-60, beta 61-66
    format "%-6s%5s %-4s %-4sA%4d    %8.3f%8.3f%8.3f%6s%6s      SEG" \
        $rec $serial $name $resname $resid 0.0 0.0 0.0 $occ $beta
}

# 1. A label PDB with what a program writing one may produce:
#    - group atoms on ATOM and on HETATM records,
#    - "no group" written as 0.00, as 0.0, and as -0.00,
#    - a serial above 99999 in VMD's hexadecimal form (186a1 = 100001),
#    - non-atom records (REMARK, TER, END) to skip.
set label [join [list \
    "REMARK  label PDB written for a test" \
    [pdbline ATOM 1 N ALA 1 " 14.01" "  1.00"] \
    [pdbline ATOM 2 CA ALA 1 " 12.01" "  1.00"] \
    [pdbline ATOM 3 C ALA 1 " 12.01" "  0.00"] \
    [pdbline ATOM 4 O ALA 1 " 16.00" "   0.0"] \
    [pdbline ATOM 5 CB ALA 1 " 12.01" " -0.00"] \
    [pdbline HETATM 6 C1 LIG 2 " 12.01" "  2.00"] \
    [pdbline HETATM 7 O1 LIG 2 " 16.00" "  2.00"] \
    "TER" \
    [pdbline ATOM 186a1 OH2 TIP3 3 " 16.00" "  3.00"] \
    "END"] \n]\n
set f [writefile label.pdb $label]
set serArray {}; set masses {}; set atnames {}; set na 0
set n [read_centersPDB $f serArray masses atnames na]
check "label: number of groups" $n 3
check "label: atoms read (ATOM and HETATM)" $na 8
check "label: group 1 (beta 1) members" [lindex $serArray 0] {1 2}
check "label: group 2 (beta 2, HETATM) members" [lindex $serArray 1] {6 7}
check "label: group 3 (hex serial) members" [lindex $serArray 2] {100001}
checknum "label: group 1 mass" [lindex $masses 0] 26.02 1e-9
checknum "label: group 2 mass" [lindex $masses 1] 28.01 1e-9
checknum "label: group 3 mass" [lindex $masses 2] 16.00 1e-9

# 2. Label PDBs that must be refused, with a message, not misread.
proc refuses {label text pattern} {
    set f [writefile bad.pdb $text]
    set serArray {}; set masses {}; set atnames {}; set na 0
    set rc [catch {read_centersPDB $f serArray masses atnames na} msg]
    check "$label: refused" $rc 1
    check "$label: message" [string match $pattern $msg] 1
}
refuses "gap in group numbers" [join [list \
    [pdbline ATOM 1 N ALA 1 " 14.01" "  1.00"] \
    [pdbline ATOM 2 CA ALA 1 " 12.01" "  3.00"]] \n]\n "*skip group 2*"
refuses "fractional beta" [join [list \
    [pdbline ATOM 1 N ALA 1 " 14.01" "  1.50"]] \n]\n "*not a whole number*"
refuses "negative beta" [join [list \
    [pdbline ATOM 1 N ALA 1 " 14.01" " -1.00"]] \n]\n "*negative*"
refuses "missing mass" [join [list \
    [pdbline ATOM 1 N ALA 1 "      " "  1.00"]] \n]\n "*occupancy*"

# 3. cv.inp: TYPE then 0-based group indices; # comments and blank lines.
set f [writefile cv.inp "# phi and psi, then a distance\nDIHED 0 1 2 3\nDIHED 1 2 3 4\n\nBOND 0 4\n"]
set cvL {}
read_cvs $f cvL 5
check "cv.inp: number of CVs" [llength $cvL] 3
check "cv.inp: first CV" [lrange [lindex $cvL 0] 0 1] {DIHED {0 1 2 3}}
check "cv.inp: third CV" [lrange [lindex $cvL 2] 0 1] {BOND {0 4}}

# 4. restr.inp: one coefficient per CV, then {key value} options, long or
#    short key.
set f [writefile restr.inp "# one restraint per CV\n1 0 0 {k 100.0} {rf PERIODIC} {b 1.0261928232448847}\n0 1 0 {SpringConstant 50} {TargetValue -0.5}\n0 0 1 {k 10} {kT 3.0} {gamma 500} {dt 0.002}\n"]
set rL {}
set nR [read_restraints $f 3 rL]
check "restr.inp: number of restraints" $nR 3
check "restr.inp: coefficients" [restr_getcvc [lindex $rL 1]] {0 1 0}
check "restr.inp: short key k" [restr_getopt [lindex $rL 0] SpringConstant k 0] 100.0
check "restr.inp: long key" [restr_getopt [lindex $rL 1] SpringConstant k 0] 50
check "restr.inp: rf" [restr_getopt [lindex $rL 0] RestraintFunction rf HARMONIC] PERIODIC
check "restr.inp: target at full precision" [restr_getopt [lindex $rL 0] TargetValue b 0] 1.0261928232448847
check "restr.inp: TAMD kT" [restr_getopt [lindex $rL 2] TAMDkT kT -1] 3.0
check "restr.inp: TAMD gamma" [restr_getopt [lindex $rL 2] TAMDgamma gamma -1] 500
check "restr.inp: absent option takes its default" [restr_getopt [lindex $rL 1] RestraintFunction rf HARMONIC] HARMONIC

file delete -force $tmp
puts "input_files: $npass passed, $nfail failed"
exit [expr {$nfail > 0}]
