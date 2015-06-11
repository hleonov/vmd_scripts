# This script is designed to search for "proton wires" starting from a certain point
# The mechanism of the script is to search for an outer layer, within some criteria, and mark all the atoms.
# The next layer follows the same criteria but will disregard atoms already found in the previous layer.
# Thus, each step an outer layer if found, and is added to the cummulative, already found list..
# vmdt -e ~/Research/scripts/water_wire.tcl -args md_0_8.pdb 1449


mol load pdb [lindex $argv 0]
set cutoff 3.5; # cutoff in angstrems
set start_res [lindex $argv 1]
set start_atoms [atomselect top "protein and resid $start_res and sidechain and type \"O.*\" \"N.*\""]

set current [$start_atoms list];
set cummulative $current
set i 0


while {$i < [lindex $argv 2]} {
  #set front [atomselect top "(type \"O.*\" \"N.*\" and within $cutoff of (index $current)) and not (index $cummulative $current)"]    
  set front [atomselect top "(type \"O.*\" and within $cutoff of (index $current)) and not (index $cummulative $current)"]    
  set cummulative [concat $cummulative $current];
  set current [$front list];
  incr i;
  #puts "Cycle $i: ([$front num]) $current ===================";
  # check for end
  if {$current == 0} {
    break;
  }
}

puts "index $cummulative"

exit;





		


	
