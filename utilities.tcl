package require Orient

namespace import Orient::orient

# =================== General useful stuff ======================

proc dihedral {a1 a2 a3 a4 fr} {
  set sel1 [atomselect top "index $a1" frame $fr]
  set sel2 [atomselect top "index $a2" frame $fr]
  set sel3 [atomselect top "index $a3" frame $fr]
  set sel4 [atomselect top "index $a4" frame $fr]
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

proc pca {seltxt fr} {
	set sel [atomselect top "$seltxt" frame $fr]
	set I [draw principalaxes $sel]
	namespace forget Orient::orient
	return $I
}

proc center_by_selection {txtsel} {
	set sel [atomselect top "$txtsel"]
	for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
		molinfo top set frame $i; 
		$sel moveby [vecinvert [measure center $sel]];
	}		
}

proc center_frame_by_selection {txtsel fr} {
  set selection [atomselect top "$txtsel" frame $fr]
  set all [atomselect top "all" frame $fr]
  $all moveby [vecinvert [measure center $selection weight mass]]
}

proc align_trj_by_prot {fr} {
	set all [atomselect top "all"]
	set sel [atomselect top "protein and type CA"]
	set ref [atomselect top "protein and type CA" frame $fr]
	for {set frame 0} {$frame < [molinfo $mol get numframes]} {incr frame} {
      $sel frame $frame
      $all frame $frame
      $all move [measure fit $sel $ref]
    }
}

proc plane_normal {seltxt fr} {
	set I [pca $seltxt $fr]
	puts "I: $I"
	set sortI [lsort -command sort_by_z $I]
	set PCA [lreplace $sortI 2 2]
	set n [norm_to_plane [lindex $PCA 0] [lindex $PCA 1]]
	puts "normal: $n"
	return $n
}
# =============== Vectors and numerical procedures ================

proc sort_by_z {a b} {
	set az [expr abs([lindex $a 2])]
	set bz [expr abs([lindex $b 2])]
	if {$az < $bz} {
		return -1
	} elseif {$az > $bz} {
		return 1
	}
	return 0
}
proc norm_to_plane {v1 v2} {
	set cross [vecnorm [veccross $v2 $v1]]
	return $cross
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

proc rad_to_deg {angle} {
  return [expr $angle * 180 / 3.14159265]
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



 
# =========== input / output procedures ==============
proc report {frame number} {
  if { $frame % $number == 0} {
    puts "Analyzing frame: $frame\n"
  }
}

#proc read_list {input} {
#	
#	set parts [split $input "-"];
#
#	global atoms;
#	set i 0;
#	set file_id [open $list_file r];
#	#gets with 2 arguments reads into $line, and returns length, -1 for EOF
#	while {[gets $file_id line] >= 0} {
#		#split line into a list that can be accessed with lindex
#		set atoms($i) [split $line " "]; 
#		incr i;		
#	}
#
#}
proc out_file {file} {
  while { [file exists $file] } {
  #	cp $file "backup.$file"
   # puts stderr "File $file exists, overwrite \[y\]"
#    if {[gets stdin] == "y"} {
#      break
#    } else {
#      puts stderr "Please suplly output file name.\n";
#      set file [gets stdin]
    }
  #}
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












