if {![info exists CFACV_BASEDIR]} {
  if {[info exists env(CFACV_BASEDIR)]} {
    set CFACV_BASEDIR $env(CFACV_BASEDIR)
  } else {
    set HOME $env(HOME)
    set CFACV_BASEDIR ${HOME}/cfacv
  }
}

load ${CFACV_BASEDIR}/lib/libgenericdataspace.so libgenericdataspace
###########################################################
# Primitive data handling procedures
###########################################################

# ListToArray: allocates a new array and assigns it values
# from the Tcl list $l; returns pointer to new array.  List
# elements are treated as double-precision by default.
#
proc ListToArray {l} {
    set length [llength $l]
    set a [new_array $length]
    set i 0
    foreach item $l {
        array_setitem $a $i $item
        incr i 1
    }
    return $a
}

# ListToArray_Data:  Assigns elements of an existing 
# array pointed to by $a from the tcl list $l
#
proc ListToArray_Data { a l } {
    set i 0
    foreach item $l {
	array_setitem $a $i $item
	incr i 1
    }
    return $a
}

# intListToArray: like listToArray except elements of
# $l are treated as integers
#
proc intListToArray {l} {
    set length [llength $l]
    set a [new_arrayint $length]
    set i 0
    foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
    return $a
}

# intListToArray_Data: list listToArray_Data except elements
# of $l are treated as integers
#
proc intListToArray_Data {a l} {
   set length [llength $l]
   set i 0
   foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
  return $a
}

# ArrayToList:  Creates a new tcl list of length $n and
# assigns it elements of existing array $a
#
proc ArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [array_getitem $a $i]
    }
    return $l
}

# intArrayToList: like ArrayToList but treats
# elements of $a as integers
#
proc intArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [arrayint_getitem $a $i]
    }
    return $l
}

proc Tcl_GenDataSpace_Write { gds i dlist } {
    set addr [GenDataSpace_getAddr $gds $i]
    ListToArray_Data $addr $dlist
    return $addr
}

proc Tcl_IGenDataSpace_Write { ids i ilist } {
    set addr [IGenDataSpace_getAddr $ids $i]
    intListToArray_Data $addr $ilist
    return $addr
}

proc Tcl_GenDataSpace_Read { gds i N } {
    set addr [GenDataSpace_getAddr $gds $i]
    return [ArrayToList $addr $N]
}

proc Tcl_IGenDataSpace_Read { ids i N } {
    set addr [IGenDataSpace_getAddr $ids $i]
    return [intArrayToList $addr $N]
}

proc aa_321 { aa3 } {
    set aa1 {}
    foreach aa $aa3 {
        switch $aa {
            ALA {
                set a1 A
            }
            VAL {
                set a1 V
            }
            ILE {
                set a1 I
            }
            LEU {
                set a1 L
            }
            PRO {
                set a1 P
            }
            PHE {
                set a1 F
            }
            TRP {
                set a1 W
            }
            MET {
                set a1 M
            }
            GLY {
                set a1 G
            }
            SER {
                set a1 S
            }
            THR {
                set a1 T
            }
            CYS {
                set a1 C
            }
            ASN {
                set a1 N
            }
            GLN {
                set a1 Q
            }
            TYR {
                set a1 Y
            }
            HIS -
            HSD -
            HSE -
            HSP {
                set a1 H
            }
            LYS {
                set a1 K
            }
            ARG {
                set a1 R
            }
            ASP {
                set a1 D
            }
            GLU {
                set a1 E
            }
            default {
		puts "aa? $aa"
                set a1 "?"
            }
        }
        lappend aa1 $a1
    }
    return $aa1
}
