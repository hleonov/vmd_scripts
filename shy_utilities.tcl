########## bigdcd and waitfor need to be first. ##############
proc bigdcd { trj script args } {
  global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame

  set bigdcd_frame 0
  set bigdcd_firstframe [molinfo top get numframes]
  set bigdcd_proc $script

  uplevel #0 trace variable vmd_frame w bigdcd_callback
  foreach dcd $args {
    animate read $trj $dcd waitfor 0
  }
}
proc bigdcd_callback { name1 name2 op } {
  global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame

  # If we're out of frames, we're also done
  set thisframe $vmd_frame($name2)
  if { $thisframe < $bigdcd_firstframe } {
    bigdcd_done
    return
  }

  incr bigdcd_frame
  if { [catch {uplevel #0 $bigdcd_proc $bigdcd_frame} msg] } {
    puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
    bigdcd_done
    return
  }
  animate delete beg $thisframe end $thisframe
  return $msg
}
proc bigdcd_done { } {
  puts "bigdcd_done"
  after idle uplevel #0 trace vdelete vmd_frame w bigdcd_callback
}
proc waitfor_callback { name1 name2 op } {
  global waitfor_file
  global waitfor_proc

  upvar $name1 arr
  set fname $arr($name2)
  if { ! [string match $fname $waitfor_file ] } { return }

  eval uplevel #0 $waitfor_proc
  after idle uplevel #0 trace vdelete vmd_trajectory_read w waitfor_callback
}
proc waitfor { file script } {
  global waitfor_file
  global waitfor_proc

  set waitfor_file $file
  set waitfor_proc $script

  uplevel #0 trace variable vmd_trajectory_read w waitfor_callback
} 
################### BigDCD done! ##################################



proc report {frame number args} {
  if { $frame % $number == 0} {
    if {[llength $args] > 0} {
      puts "Analyzing frame: $frame [lindex $args 0]"

    } else {
      puts "Analyzing frame: $frame"
    }
  }
}
proc center_by_selection {sel} {
  set selection [atomselect top $sel]
  set all [atomselect top "all"]
  $all moveby [vecinvert [measure center $selection weight mass]]
}
proc out_file {file} {
  # if there is a command line option called "-rf"
  # then we don't check for overwriting
  global argv
  set overwrite_check "yes"
  for {set i 0} {$i < [llength $argv]} {incr i} {
    if {[lindex $argv $i] == "-rf"} {
      set overwrite_check "no"
      if { [file exists $file] } {
        puts stderr "Overwriting $file\n";
      }
    }
  }
  if {$overwrite_check == "yes"} {
    while { [file exists $file] } {
      puts stderr "File $file exists, overwrite \[y\]"
      if {[gets stdin] == "y"} {
    	break
      } else {
    	puts stderr "Please suplly output file name.\n";
    	set file [gets stdin]
      }
    }
  }
  return $file
}
proc check_for_help {message} {
  global argv
  for {set i 0} {$i < [llength $argv]} {incr i} {
    set x [lindex $argv $i]
    if {$x == "-h" || $x == "-help" || $x == "-H" || $x == "-Help"} {
      puts stderr "\n\n$message\n\n"
      exit
    }
  }
}
proc angle_between_two_vectors { a b } {
  # Angle between two vectors
  set amag [veclength $a]
  set bmag [veclength $b]
  set dotprod [vecdot $a $b]
  set normalized [expr {$dotprod / ($amag * $bmag)}]
  if {$normalized < -1} {
    set normalized -1
  }
  if {$normalized > 1} {
    set normalized 1
  }
  return [expr {57.2958 * acos($normalized)}]
}
proc max {data} {
  return [lindex [lsort $data] end]
}
proc min {data} {
  return [lindex [lsort $data] 0]
}
proc histogram {data bin_number min max} {
  # this will generate a histogram according to the
  # supplied number of bin of the given data with the given min max
  # it will retrun it as 2 lists: x and y
  set range    [expr $max - $min]
  set bin_size [expr 1.0 * $range / $bin_number]
  set x {}
  set y {}
  # construct the x axis
  for {set i $min} {$i <= $max} {set i [expr $i + $bin_size]} {
    lappend x $i
    lappend y 0
  }
  # decide which bin is it in
  foreach element $data {
    for {set i 0} {$i <= $bin_number} {incr i} {
      if {$element >= [lindex $x $i] && $element < [lindex $x [expr $i + 1]]} {
	lset y $i [expr [lindex $y $i] + 1] 
      }
    }
  }
  return [list $x $y]
}
proc average_list {list_a list_b ratio} {
  # this will average two lists (of the same length)
  # using a ratio. Thus if the ratio is > 1
  # then list_a will be over representative
  # It will retrun the averaged list
  set average {}
  for {set i 0} {$i < [llength $list_a]} {incr i} {
    set x [expr ($ratio * [lindex $list_a $i] + [lindex $list_b $i]) / (2.0 + ($ratio - 1))];
    lappend average $x;
  }
  return $average
}
proc add_number_to_list {input_list number} {
  for {set i 0} {$i < [llength $input_list]} {incr i} {
    lset input_list $i [expr [lindex $input_list $i] + $number]
  }
  return $input_list
}
proc vmd_draw_arrow {mol start end thickness_scaling} {
    set scaling [expr [veclength [vecsub $end $start]]/100]
    set scaling [expr $scaling * $thickness_scaling]
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius [expr 2*$scaling]
    puts [list cone $middle $end radius [expr 5*$scaling]]
    graphics $mol cone $middle $end radius [expr 5*$scaling]
}
proc rotation_x {angle} {
  # rotate a molecule $angle degress (in deg) about the x axis
  set all [atomselect top "all"];
  set matrix [transaxis x $angle deg];
  $all move $matrix
  #return $matrix
}
proc rotation_y {angle} {
  # rotate a molecule $angle degress (in deg) about the y axis
  set all [atomselect top "all"];
  set matrix [transaxis y $angle deg];
  $all move $matrix
  #return $matrix
}

proc rotation_z {angle} {
  # rotate a molecule $angle degress (in deg) about the z axis
  set all [atomselect top "all"];
  set matrix [transaxis z $angle deg];
  $all move $matrix
  #return $matrix
}




































