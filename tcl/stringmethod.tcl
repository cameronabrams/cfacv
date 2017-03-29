# String method using NAMD's built-in replica implementation
# and the CFACV tclforces module
#
# Cameron F Abrams
# 2014-2016

# check for base directory name environment variable;
# if not set, use default
if {![info exists CFACV_BASEDIR]} {
  if {[info exists env(CFACV_BASEDIR)]} {
    set CFACV_BASEDIR $env(CFACV_BASEDIR)
  } else {
    set HOME $env(HOME)
    set CFACV_BASEDIR ${HOME}/research/cfacv
  }
}

# will die if CFACV_BASEDIR is invalid
source ${CFACV_BASEDIR}/tcl/cfacv.tcl

replicaBarrier

set nr [numReplicas]
set r [myReplica]
set replica_id $r

puts stdout "CFACV/SM) replica id $replica_id"

# check for restart
# $restart_root is set by the restart file
# sourced in jobN.conf (N>0) 
if {[info exists restart_root]} {
  set restart_root [format $restart_root $replica_id]
  source $restart_root.$replica_id.tcl
} else {
  set i_job 0
  set i_run 0
  set i_step 0
  if {[info exists first_timestep]} {
    set i_step $first_timestep
  }

  set replica(index) $r
  set replica(is_Reactant) [expr {$r == 0}]
  set replica(is_Product) [expr {($r + 1) == $nr}]
  if {!$replica(is_Reactant)} {
      set replica(back) [expr $r - 1]
  } else {
      set replica(back) NULL
  }
  if {!$replica(is_Product)} {
      set replica(forward) [expr $r + 1]
  } else {
      set replica(forward) NULL
  } 
}

set job_output_root "$output_root.job$i_job"
set job_initial_root "$initial_root.job$i_job"

# if i_run > 0, this is a restart
if {$i_run} {
  bincoordinates $restart_root.$replica_id.coor
  binvelocities $restart_root.$replica_id.vel
  extendedSystem $restart_root.$replica_id.xsc
} else {
  bincoordinates $job_initial_root.$replica_id.coor
  if { [file exists $job_initial_root.$replica_id.vel] } {
     binvelocities $job_initial_root.$replica_id.vel
  } else {
     temperature $TEMP
  }
  extendedSystem $job_initial_root.$replica_id.xsc
}

firsttimestep $i_step

langevinTemp $TEMP
outputname [format $job_output_root.$replica_id $replica_id]
outputEnergies $outputEnergiesEvery
dcdFreq [expr $steps_per_run * $runs_per_frame]

set SMPARAMS(nu) 0.0

source $namd_config_file

set history_file [open [format "$job_output_root.$replica_id.history" $replica_id] "w"]
fconfigure $history_file -buffering line

set SMPARAMS(evolve) 1
if { $replica(is_Reactant) || $replica(is_Product) } {
    if {!$SMPARAMS(evolve_ends)} {
	set SMPARAMS(evolve) 0
    }
}

# if SMPARAMS(nu) is set, then enable climbing
set SMPARAMS(climb) 0
if { [expr $SMPARAMS(nu) > 0.0] } {
   set SMPARAMS(climb) 1
   puts "CFACV/SM) Climbing is enabled."
}

if { $SMPARAMS(climb) && $replica(is_Reactant) } {
   set SMPARAMS(stepsize) 0 ; # don't let the reactant image evolve
}


# read in the z-values for all images
# from the restr.inp files used to make them
# THIS MUST BE FIXED FOR RESTART!
if { [expr $replica_id == 0] } {
  global z
  set image_z(0) {}
  for {set i 0} {$i < $nr} {incr i} {
    set image_z($i) {}
    set i_restrfn [format $init_restr_file $i]
    set inStream [open $i_restrfn "r"]
    set data [read -nonewline $inStream]
    close $inStream
    set lines [split $data \n]
    foreach line $lines {
      foreach elem $line {
        if {[llength $elem] > 1} {
          set key [lindex $elem 0]
          set val [lindex $elem 1]
          if {$key == "b"||$key == "z"} {
            lappend image_z($i) $val
          }
        }
      }
    }
    puts stdout "CFACV/SM) image $i z(0) = $image_z($i) from $i_restrfn"
  }
  set sm 0
  set z $image_z(0)
  set SMPARAMS(nz) [llength $image_z(0)]
  puts stdout "CFACV/SM) CV space is $SMPARAMS(nz) dimensional"
  for {set i 1} {$i < $nr} {incr i} {
    replicaSend $image_z($i) $i
  }
} else {
  global z
  set z [replicaRecv 0]
  set SMPARAMS(nz) [llength $z]
  puts stdout "CFACV/SM) replica $replica_id z = $z dim $SMPARAMS(nz)"
}

# ALLOCATE STRING-METHOD DATASPACE on MASTER ONLY
# (all replicas have a restraint dataspace too)
if {[expr $replica_id == 0]} {
  puts stdout "CFACV/SM) Calling for SM dataspace allocation for $nr images in $SMPARAMS(nz)-dim CV SPACE"
  set smds [New_stringMethod_Dataspace $nr $SMPARAMS(nz) $SMPARAMS(outputlevel) $SMPARAMS(nu) $SMPARAMS(evolve_ends)]
  for {set i 0} {$i < $nr} {incr i} {
    ListToArray_Data [SMDataSpace_image_z $smds $i] $image_z($i)
  }
  SMDataSpace_set_reparam_tol $smds $SMPARAMS(reparam_tol) $SMPARAMS(maxiter)
}

set new_z 1

while {$i_run < $num_runs} {

  puts stdout "CFACV/SM) run $i_run begins"
  flush stdout

  # run MD steps to accumulate dF/dz and M on this replica
  run $steps_per_run

  # acquire this replica's dF/dz and metric tensor as TcL lists
  DataSpace_Tally $ds
  set g [ArrayToList [DataSpace_g $ds] $SMPARAMS(nz)]
  set M [ArrayToList [DataSpace_MM $ds] [expr $SMPARAMS(nz)*$SMPARAMS(nz)]]

  # master collects dF/dz and M from each replica
  if { [expr $replica_id == 0] } {
    set image_g(0) $g
    set image_M(0) $M
    for {set i 1} {$i < $nr} {incr i} {
       set image_g($i) [replicaRecv $i]
       set image_M($i) [replicaRecv $i]
    }
  } else {
    replicaSend $g 0
    replicaSend $M 0
  }
  
  replicaBarrier
 
  # master transfers dF/dz and M from all replicas to SM dataspace
  if {[expr $replica_id == 0]} {
    for {set i 0} {$i < $nr} {incr i} { 
       ListToArray_Data [SMDataSpace_image_g $smds $i] $image_g($i)
       ListToArray_Data [SMDataSpace_image_M $smds $i] $image_M($i)
       if { [expr ($i_run%$runs_per_frame) == 0] } {
         puts $history_file "SM iter $i_run replica $i : gradient dF/dz $image_g($i)"
         puts $history_file "SM iter $i_run replica $i : metric tensor $image_M($i)"
       }
    }

    if { $SMPARAMS(evolve) } {
      if { $SMPARAMS(climb) } {
          # if requested, master alters the gradient at the climbing end of the string
          SMDataSpace_climb $smds
      }
      # master moves all image z's
      SMDataSpace_MoveString $smds $SMPARAMS(stepsize)
    }

    # master retrieves the new z-positions from the SM dataspace and sends them to all replicas
    for {set i 0} {$i < $nr} {incr i} {
      set image_z($i) [ArrayToList [SMDataSpace_image_z $smds $i] $SMPARAMS(nz)]
      if { [expr $i > 0] } {
        replicaSend $image_z($i) $i
      } else {
         set z $image_z($i)
         set new_z 1
      }
    }
  } else {
    # this non-master replica receives the new z values
    set z [replicaRecv 0]
    set new_z 1
  }

  #set z_val $z
  # This replica resets its dF/dz and M accumulators
  Tcl_DataSpace_resetAccumulators $ds
  
  replicaBarrier
#----------------------
                          
  incr i_step $steps_per_run

  if {[expr $replica_id == 0]} {
    for {set i 0} {$i < $nr} {incr i} {
      if { [expr ($i_run%$runs_per_frame) == 0] } {
	puts $history_file "SM iter $i_run replica $i : reparam z $image_z($i)"
      }
    } 
  }

  incr i_run

  #restart
  if { $i_run % ($runs_per_frame * $frames_per_restart) == 0 ||
        $i_run == $num_runs } {  # restart
    set restart_root "$job_output_root.restart$i_run"
    output [format $restart_root.$replica_id $replica_id]
    set rfile [open [format "$restart_root.$replica_id.tcl" $replica_id] "w"]
    puts $rfile [list array set replica [array get replica]]
    close $rfile
    replicaBarrier
    if { $replica_id == 0 } {
      set rfile [open [format "$restart_root.tcl" ""] "w"]
      puts $rfile [list set i_job [expr $i_job + 1]]
      puts $rfile [list set i_run $i_run]
      puts $rfile [list set i_step $i_step]
      puts $rfile [list set restart_root $restart_root]
      close $rfile
      if [info exists old_restart_root] {
        set oldroot [format $old_restart_root ""]
        file delete $oldroot.tcl
      }
    }
    replicaBarrier
    if [info exists old_restart_root] {
      set oldroot [format $old_restart_root $replica_id]
      file delete $oldroot.$replica_id.tcl
      file delete $oldroot.$replica_id.coor
      file delete $oldroot.$replica_id.vel
      file delete $oldroot.$replica_id.xsc
    }
    set old_restart_root $restart_root
  }


}

replicaBarrier 
# output final iteration                                                                                                                
if {[expr $replica_id == 0]} {                                                                                                         
  for {set i 0} {$i < $nr} {incr i} {                                                                                                
     puts $history_file "SM FINALiter $i_run replica $i : reparam z $image_z($i)"                                                    
  }                                                                                                                                  
}

replicaBarrier

