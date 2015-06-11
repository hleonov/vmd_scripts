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
        puts stderr "Overwriting $file\n";
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
proc average_between_points {vec1 vec2} {
  set x1 [lindex $vec1 0]
  set y1 [lindex $vec1 1]
  set z1 [lindex $vec1 2]
  set x2 [lindex $vec2 0]
  set y2 [lindex $vec2 1]
  set z2 [lindex $vec2 2]
  set x [expr ($x1 + $x2) / 2]
  set y [expr ($y1 + $y2) / 2]
  set z [expr ($z1 + $z2) / 2]
  return [list $x $y $z]
}

proc rotate_helix {angle number_of_amino_acids_to_average buffer_from_end} {
  # this will roatte a helix about its axis by $angle degrees (in deg)
  # the helix axis is determined by averaging a number of Ca atoms
  # given by $number_of_amino_acids_to_average
  # it wont take into acount $buffer_from_end residues from the end
  set all [atomselect top all];
  set CAs [atomselect top "name CA"];
  # find the range of residues 
  set first_resid [lindex [$CAs get {resid}] 0];
  set last_resid  [lindex [$CAs get {resid}] end];
  # calculate the helix axis
  set a [expr $first_resid + $buffer_from_end]
  set b [expr $first_resid + $buffer_from_end + $number_of_amino_acids_to_average]
  set start [atomselect top "name CA and resid $a to $b"];
  set a [expr $last_resid - $buffer_from_end - $number_of_amino_acids_to_average]
  set b [expr $last_resid - $buffer_from_end]
  set end  [atomselect top "name CA and resid $a to $b"];
  set start_coords [measure center $start];
  set end_coords   [measure center $end  ];
  set matrix [trans bond $start_coords $end_coords $angle deg]
  $all move $matrix
}

proc rotate_x {angle} {
  # rotate a molecule $angel degress (in deg) about the x axis
  set all [atomselect top "all"];
  set matrix [transaxis x $angle deg];
  $all move $matrix
}
proc rotate_y {angle} {
  # rotate a molecule $angel degress (in deg) about the y axis
  set all [atomselect top "all"];
  set matrix [transaxis y $angle deg];
  $all move $matrix
}
proc rotate_z {angle} {
  # rotate a molecule $angel degress (in deg) about the z axis
  set all [atomselect top "all"];
  set matrix [transaxis z $angle deg];
  $all move $matrix
}

proc rotation_x {angle} {
  # rotate a molecule $angel degress (in deg) about the x axis
  set all [atomselect top "all"];
  set matrix [transaxis x $angle deg];
  # $all move $matrix
  return $matrix
}
proc rotation_y {angle} {
  # rotate a molecule $angel degress (in deg) about the y axis
  set all [atomselect top "all"];
  set matrix [transaxis y $angle deg];
  # $all move $matrix
  return $matrix
}

proc rotation_z {angle} {
  # rotate a molecule $angel degress (in deg) about the z axis
  set all [atomselect top "all"];
  set matrix [transaxis z $angle deg];
  # $all move $matrix
  return $matrix
}

proc rotation_helix {angle number_of_amino_acids_to_average buffer_from_end} {
  # this will roatte a helix about its axis by $angle degrees (in deg)
  # the helix axis is determined by averaging a number of Ca atoms
  # given by $number_of_amino_acids_to_average
  # it wont take into acount $buffer_from_end residues from the end
  set all [atomselect top all];
  set CAs [atomselect top "name CA"];
  # find the range of residues 
  set first_resid [lindex [$CAs get {resid}] 0];
  set last_resid  [lindex [$CAs get {resid}] end];
  # calculate the helix axis
  set a [expr $first_resid + $buffer_from_end]
  set b [expr $first_resid + $buffer_from_end + $number_of_amino_acids_to_average]
  set start [atomselect top "name CA and resid $a to $b"];
  set a [expr $last_resid - $buffer_from_end - $number_of_amino_acids_to_average]
  set b [expr $last_resid - $buffer_from_end]
  set end  [atomselect top "name CA and resid $a to $b"];
  set start_coords [measure center $start];
  set end_coords   [measure center $end  ];
  set matrix [trans bond $start_coords $end_coords $angle deg]
  return $matrix
  # $all move $matrix
}


proc rad_to_deg {angle} {
  return [expr $angle * 180 / 3.14159265358979323846]
}
proc dihedral {a1 a2 a3 a4} {
  set sel1 [atomselect top "index $a1"]
  set sel2 [atomselect top "index $a2"]
  set sel3 [atomselect top "index $a3"]
  set sel4 [atomselect top "index $a4"]
  set coord1 [lindex [$sel1 get {x y z}] 0]
  set coord2 [lindex [$sel2 get {x y z}] 0]
  set coord3 [lindex [$sel3 get {x y z}] 0]
  set coord4 [lindex [$sel4 get {x y z}] 0]
  set v1 [vecsub $coord1 $coord2] 
  set v2 [vecsub $coord3 $coord2]
  set v3 [vecsub $coord4 $coord3]
  set cross1 [vecnorm [veccross $v2 $v1]]
  set cross2 [vecnorm [veccross $v2 $v3]]
  set dot [vecdot $cross1 $cross2]
  set angle [expr acos($dot)]
  return [rad_to_deg $angle]
}

proc dihedral_z {a1 a2 a3} {
  set sel1 [atomselect top "index $a1"]
  set sel2 [atomselect top "index $a2"]
  set sel3 [atomselect top "index $a3"]
  set coord1 [lindex [$sel1 get {x y z}] 0]
  set coord2 [lindex [$sel2 get {x y z}] 0]
  set coord3 [lindex [$sel3 get {x y z}] 0]
  set coord4 {0 0 1};
  set v1 [vecsub $coord1 $coord2] 
  set v2 [vecsub $coord3 $coord2]
  set v3 [vecsub $coord4 $coord3]
  set cross1 [vecnorm [veccross $v2 $v1]]
  set cross2 [vecnorm [veccross $v2 $v3]]
  set dot [vecdot $cross1 $cross2]
  set angle [expr acos($dot)]
  return [rad_to_deg $angle]
}




















