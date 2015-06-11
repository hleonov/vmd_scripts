# this script should run through a trajectory and generate a histogram with waters
# the outcome is waters per A cubed. Bulck water density is 0.033.

# it runs on a single xtc and takes a while since it uses bigdcd to read one frame at a time

# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/analyze_waters_all.tcl -args system_from_em.pdb  system.xtc results_163.xls results_164.xls 140 35.0 3 155 156



source ~arkini/data/scripts/Tcl/bigdcd.tcl
source ~arkini/data/scripts/Tcl/waitfor.tcl

# the command line arguments parsing
set pdb               [lindex $argv 0] ;# system_from_em_free.pdb
set xtc               [lindex $argv 1] ;# system.xtc
set file_a            [lindex $argv 2] ;# results_163.xls
set file_b            [lindex $argv 3] ;# results_164.xls
set number_of_slices  [lindex $argv 4] ;# 50
set range             [lindex $argv 5] ;# 35
set radius            [lindex $argv 6] ;# 3
set residue_a         [lindex $argv 7] ;# 155
set residue_b         [lindex $argv 8] ;# 156

# here I load the protein  
mol load pdb $pdb  
# define some basic parameters
set framestotal 0
set total_number_frames 0
# some simple calculations
set slab_width [expr 2 * $range * 1.0 / $number_of_slices ]
set r2 [expr $radius*$radius]
# open the file handle
set fid_a [open $file_a w]
set fid_b [open $file_b w]
# initialze the counting array
for {set i 0} {$i <= $number_of_slices + 1} {incr i} {
  set counting_array_a($i) 0
  set counting_array_b($i) 0
}
# run thru the trajectory
proc analyze {frame} {
	
  global counting_array_a counting_array_b r2 range slab_width total_number_frames residue_a residue_b
  
  # the following centers the system aroudn the selection and then shifts is by the vertical range
  # this is just for ease of computation. # note that the cylinder is now centered at the selection at EACH frame.
  # for the 1st residue
  set all [atomselect top "all"]
  set sel [atomselect top "resid $residue_a and protein and type OD1 OD2"]
  $all moveby [vecinvert [measure center $sel]]
  $all moveby [list 0 0 $range]
  set water_sel [atomselect top "name OWS and x^2 + y^2 <= $r2 and (z >= 0 ) and (z <= 2* $range )"]
  # now count the number of atoms in each histogram bin
  foreach z [$water_sel get {z}] {     
    set slice [ expr floor( [expr $z / $slab_width ] ) + 1 ] 
    set slice [ expr int($slice)] 
    incr counting_array_a($slice)
  }
  # for the 2nd residue
  set sel [atomselect top "resid $residue_b and protein and type OD1 OD2"]
  $all moveby [vecinvert [measure center $sel]]
  $all moveby [list 0 0 $range]
  set water_sel [atomselect top "name OWS and x^2 + y^2 <= $r2 and (z >= 0 ) and (z <= 2* $range )"]
  # now count the number of atoms in each histogram bin
  foreach z [$water_sel get {z}] {     
    set slice [ expr floor( [expr $z / $slab_width ] ) + 1 ] 
    set slice [ expr int($slice)] 
    incr counting_array_b($slice)
  }
  # this makes a nice printout every 25 frames
  if {[expr $total_number_frames / 25.0] == int([expr $total_number_frames / 25.0])} {
   puts "Analyzing frame: $total_number_frames"
  }
  incr total_number_frames
}

proc end {} {
  global total_number_frames r2 slab_width fid_a fid_b number_of_slices counting_array_a counting_array_b
  set slab_volume [expr $r2 * 3.14159265358979 * $slab_width]
  for {set l 1} {$l < $number_of_slices} {incr l} { 
    set counting_array_a($l) [expr $counting_array_a($l) * 1.0 / ( $total_number_frames * $slab_volume ) ]
    puts $fid_a "$l $counting_array_a($l)"
  }
  for {set l 1} {$l < $number_of_slices} {incr l} { 
    set counting_array_b($l) [expr $counting_array_b($l) * 1.0 / ( $total_number_frames * $slab_volume ) ]
    puts $fid_b "$l $counting_array_b($l)"
  }
  close $fid_a
  close $fid_b
  exit
}


waitfor $xtc end
bigdcd xtc analyze $xtc

