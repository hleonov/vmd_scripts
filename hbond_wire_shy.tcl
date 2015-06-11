# This script is designed to search for "proton wires" starting from a certain point
# The mechanism of the script is to search for an outer layer, within some criteria, and mark all the atoms.
# The next layer follows the same criteria but will disregard atoms already found in the previous layer.
# Thus, each step an outer layer if found, and is added to the cummulative, already found list..
# vmdt -e ~/Research/scripts/water_wire.tcl -args md_0_8.pdb 155 12 md_0_8.xtc results.txt 5


if {![file exists [lindex $argv 0]]} {
  puts "Please supply pdb as 1st argumnet";
  exit;
} 
if {![string is integer -strict [lindex $argv 1]]} {
  puts "Please supply resid as 2nd argumnet";
  exit;
}
if {![string is integer -strict [lindex $argv 2]]} {
  puts "Please supply number fo repititions as 3rd argumnet";
  exit;
}
if {![file exists [lindex $argv 3]]} {
  puts "Please supply xtc as 4th argumnet";
  exit;
} 
if {![lindex $argv 4]} {
  puts "Please supply name fo output file as 5th argumnet";
  exit;
} 
if {![string is integer -strict [lindex $argv 5]]} {
  puts "Please supply number frames to skip as 6th argumnet";
  exit;
}


mol load pdb [lindex $argv 0]
animate read xtc [lindex $argv 3] skip [lindex $argv 5] waitfor all
set cutoff 3.5; # cutoff in angstrems
set start_res [lindex $argv 1]
# if you want to use N
set use_n 1
set fid [open [lindex $argv 4] w]


proc get_h_bond_list {frame} {
  global start_res cutoff use_n argv
  set start_atoms [atomselect top "protein and resid $start_res and sidechain and type \"O.*\" \"N.*\"" frame $frame]
  set current [$start_atoms list];
  set cummulative $current
  set i 0
  while {$i < [lindex $argv 2]} {
    if {$use_n} {
      set front [atomselect top "(type \"O.*\" \"N.*\" and within $cutoff of (index $current)) and not (index $cummulative $current)"  frame $frame]    
    } else {
      set front [atomselect top "(type \"O.*\" and within $cutoff of (index $current)) and not (index $cummulative $current)"  frame $frame]    
    }
    set cummulative [concat $cummulative $current];
    set current [$front list];
    incr i;
    if {$current == 0} {
      break;
    }
  }
  return $cummulative;
}

for {set i 1} {$i < [molinfo top get numframes]} {incr i} {
  puts "analyzing frame $i"
  puts $fid [get_h_bond_list $i];
}





exit;





		


	
