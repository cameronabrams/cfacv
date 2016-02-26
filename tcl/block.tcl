if {![info exists CFACV_BASEDIR]} {
  if {[info exists env(CFACV_BASEDIR)]} {
    set CFACV_BASEDIR $env(CFACV_BASEDIR)
  } else {
    set HOME $env(HOME)
    set CFACV_BASEDIR ${HOME}/cfacv
  }
}

source $CFACV_BASEDIR/tcl/genericdataspace.tcl

proc block_by_seq { selList molid selstr nBlock } {
    upvar $selList p
    global ids
    global gds
    
    set chop_me [atomselect $molid "$selstr and name CA"]
    set resids [$chop_me get resid]
    set N [$chop_me num]
    set npp [expr $N / $nBlock ]
    set rr [expr $N - ($nBlock * $npp)]
    $chop_me delete
    puts "CFACV) blocking $N residues ([lindex $resids 0] to [lindex $resids [expr [llength $resids] - 1]]) into $nBlock blocks gives $rr remainder"

    set bsize {}
    set is_one_left [expr $rr > 0]
    for {set i 0} { $i < $nBlock } { incr i } {
	lappend bsize [expr $npp + $is_one_left]
	if { $is_one_left } {
	    set rr [expr $rr - 1]
	    set is_one_left [expr $rr > 0]
	}
    }

    set firstindex 0
    for { set i 0 } { $i < $nBlock } {incr i } {
	puts "resid [lindex $resids $firstindex] to [lindex $resids [expr $firstindex + [lindex $bsize $i]-1]]"
	set a [atomselect $molid "($selstr) and resid [lindex $resids $firstindex] to [lindex $resids [expr $firstindex + [lindex $bsize $i] - 1]]"]
	set firstindex  [expr $firstindex + [lindex $bsize $i]]
	$a global
	lappend p $a 
    }

    return $nBlock
}

proc block_by_rgyr { selList molid selstr nBlock } {
    upvar $selList p
    global ids
    global gds

    set first_bin [llength $p]

    set all [atomselect $molid "all"]
    set savebeta [$all get beta]
    $all set beta 0.0

    set sort_me [atomselect $molid "$selstr and name CA"]
    set N [$sort_me num]
    set npp [expr $N / $nBlock ]
    puts "CFACV) blocking $N residues into $nBlock blocks"

    # assign initial bins
    set nLeftOvers [expr $N-( $npp * $nBlock)]
    set bin 1
    set ibin 0
    set binlist {}
    set rglist {}
    set thisnpp $npp
    if {$nLeftOvers > 0} { 
	incr thisnpp 
	incr nLeftOvers -1
    }
    for { set i 0 } { $i < $N } { incr i } {
	set thebin $bin
	incr ibin
	if {$ibin == $thisnpp} {
	    incr bin
	    set ibin 0
	    set thisnpp $npp
	    if {$nLeftOvers > 0} { 
		incr thisnpp 
		incr nLeftOvers -1
	    }
	}
#	puts "DB: CA $i into bin $thebin"
	if {$thebin > $nBlock} {
	    set thebin $nBlock
	}
	lappend rglist 0.0
	lappend binlist $thebin
    }

    set binarr [Tcl_IGenDataSpace_Write $ids 0 $binlist]

    set xarr [Tcl_GenDataSpace_Write $gds 0 [$sort_me get x]]
    set yarr [Tcl_GenDataSpace_Write $gds 1 [$sort_me get y]]
    set zarr [Tcl_GenDataSpace_Write $gds 2 [$sort_me get z]]
    set rgarr [Tcl_GenDataSpace_Write $gds 3 $rglist]

    bin_sort $binarr $xarr $yarr $zarr $N $nBlock 1000000 [clock seconds]

    set binlist [Tcl_IGenDataSpace_Read $ids 0 $N]
    set rglist [Tcl_GenDataSpace_Read $gds 3 $nBlock]

#    puts "TCL: after bin_sort: bin $binlist"
#    puts "TCL: first_bin $first_bin"

    set i 0
    foreach r [$sort_me get resid] {
	set a [atomselect $molid "($selstr) and resid $r"]
	set thisBin [expr [lindex $binlist $i] + $first_bin]
#	puts "TCL: setting beta of resid [lsort -unique [$a get resname]][lsort -unique [$a get resid]] to $thisBin"
	$a set beta $thisBin
	$a delete
	incr i
    }
    
    for { set i 0 } { $i < $nBlock } {incr i } {
	set a [atomselect $molid "($selstr) and beta [expr $first_bin + $i + 1]"]
	puts "CFACV) Beta [expr $first_bin + $i + 1] owns [$a num] atom-members in [llength [lsort -unique [$a get resid]]] residues."
	$a global
	lappend p $a 
    }

    $all set beta $savebeta
    $all delete

    return $nBlock
}

proc binary_block_by_rgyr { selList molid selstr minNRes } {
    upvar $selList p
    global ids
    global gds

    set first_bin [llength $p]

    set all [atomselect $molid "all"]
    set savebeta [$all get beta]
    $all set beta 0.0

    set sort_me [atomselect $molid "$selstr and name CA"]
    set N [$sort_me num]
#    set npp [expr $N / $np ]
#    puts "$N $npp $np : [expr $npp * $np]"

    # assign initial bins
#    set nLeftOvers [expr $N-( $npp * $np)]
#    set bin 1
#    set ibin 0
    set binlist {}
    set rglist {}
#    set thisnpp $npp
#    if {$nLeftOvers > 0} { 
#	incr thisnpp 
#	incr nLeftOvers -1
#    }
    for { set i 0 } { $i < $N } { incr i } {
#	set thebin $bin
#	incr ibin
#	if {$ibin == $thisnpp} {
#	    incr bin
#	    set ibin 0
#	    set thisnpp $npp
#	    if {$nLeftOvers > 0} { 
#		incr thisnpp 
#		incr nLeftOvers -1
#	    }
#	}
#	puts "DB: CA $i into bin $thebin"
#	if {$thebin > $np} {
#	    set thebin $np
#	}
	lappend binlist 0 ; # only if using recursion
	lappend rglist 0.0
#	lappend binlist $thebin
    }

    set binarr [Tcl_IGenDataSpace_Write $ids 0 $binlist]

    set xarr [Tcl_GenDataSpace_Write $gds 0 [$sort_me get x]]
    set yarr [Tcl_GenDataSpace_Write $gds 1 [$sort_me get y]]
    set zarr [Tcl_GenDataSpace_Write $gds 2 [$sort_me get z]]
    set rgarr [Tcl_GenDataSpace_Write $gds 3 $rglist]

    #bin_sort $binarr $xarr $yarr $zarr $N $np 1000000 [clock seconds]
    set np [rgyr_sort [Null_centerStruct] $binarr $xarr $yarr $zarr $N $minNRes $rgarr [clock seconds]]

    set binlist [Tcl_IGenDataSpace_Read $ids 0 $N]
    set rglist [Tcl_GenDataSpace_Read $gds 3 $np]

#    puts "TCL: after bin_sort: bin $binlist"
#    puts "TCL: first_bin $first_bin"
    puts "CFACV) center Rg's: $rglist"

    set i 0
    foreach r [$sort_me get resid] {
	set a [atomselect $molid "($selstr) and resid $r"]
	set thisBin [expr [lindex $binlist $i] + $first_bin]
#	puts "TCL: setting beta of resid [lsort -unique [$a get resname]][lsort -unique [$a get resid]] to $thisBin"
	$a set beta $thisBin
	$a delete
	incr i
    }
    
    for { set i 0 } { $i < $np } {incr i } {
	set a [atomselect $molid "($selstr) and beta [expr $first_bin + $i + 1]"]
	puts "CFACV) Beta [expr $first_bin + $i + 1] owns [$a num] atom-members in [llength [lsort -unique [$a get resid]]] residues."
	$a global
	lappend p $a 
    }

    $all set beta $savebeta
    $all delete

    return $np
}
