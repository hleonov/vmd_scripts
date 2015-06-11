# Runs through a trajectory file (xtc or trr) and generates a histogram with water average over time. 
# The outcome is waters per A cubed. Bulk water density is 0.033.
# The script loads all the trajectory into memory, so when large trr exists - don't use it.

# In addition, an implementation difference from Shy: 
#	1. move the first frame protein (rather than entire system) to 0 0 0, adjust protein by $range on the Z-axis
#	2. Adjust subsequent frames by aligning to the first reference protein (measure fit), instead of #1 transformation.

# How to run: vmd -dispdev text -e ~/vmd_scripts/water.tcl -args (number of slices) (gro) (xtcfile) (skip) (trj-type) (radius)

source ~/vmd_scripts/shy_utilities.tcl
check_for_help "vmd -dispdev text -e ~/vmd_scripts/water.tcl -args 50 system_from_PR.gro system_from_MD.xtc 10 xtc 7"

set number_of_slices 	[lindex $argv 0];
set pdb			 		[lindex $argv 1];
set xtc 				[lindex $argv 2];
set skip_fr 			[lindex $argv 3]; 
set trj 				[lindex $argv 4];
set radius 				[lindex $argv 5];	
#controls the radius of the sphere going through the pore

set beg "xtc_results_"
set end "_density.dat"
set file_name [concat $beg$number_of_slices$end]

# here I load the protein  
mol load pdb $pdb
set all [atomselect top "all" frame 0]
set z_min [lindex [lindex [measure minmax $all] 0 ] 2]
set z_max [lindex [lindex [measure minmax $all] 1 ] 2]

#regular range is abs($z_max - $z_min)+1)/2. 
#I add +3 for future (xtc) water movement
set range [expr ((abs($z_max - $z_min)+1)/2) + 3]
  
# define some basic parameters
set pi 3.14159265358979

set total_number_frames 0
set slab_width [expr 2 * $range * 1.0 / $number_of_slices ]
set r2 [expr $radius*$radius]

# open the file handle
set fid [open $file_name w]
# initialze the counting array
for {set i 0} {$i <= $number_of_slices+1} {incr i} {
  set counting_array($i) 0
} 
# read xtc.
animate read $trj $xtc skip $skip_fr waitfor all
set number_of_frames [molinfo top get numframes]

# Center around the protein's center
#first frame
set start [atomselect top "all" frame 0]							
set protein [atomselect top "protein" frame 0]
#vector from center of protein to 0 0 0
set vecinv [vecinvert [measure center $protein]]					
#move first frame there
$start moveby [vecinvert [measure center $protein]]					

# Move it such that the minimal z coordiante (of future waters) is at 0
$start moveby [list 0 0 $range]
# set the reference for further alignment
set ref [atomselect top "protein and name CA" frame 0]
for {set k 1} {$k < $number_of_frames} {incr k} { 
  # Align and move
  set all [atomselect top "all" frame $k]
  set ca [atomselect top "protein and name CA" frame $k]
  $all move [measure fit $ca $ref]	
  set water_sel [atomselect top "name OW and (x^2 + y^2 <= $r2) and (z <= 2*$range)" frame $k]
  # now count the number of atoms in each histogram bin
  foreach z [$water_sel get {z}] {     
	set slice [ expr floor( [expr $z / $slab_width ] ) +1 ] 
	set slice [expr int($slice)]
	incr counting_array($slice)
  }
  incr total_number_frames
}

# calculate the slab volume
set slab_volume [expr $r2 * $pi * $slab_width]
puts $total_number_frames
for {set l 1} {$l <= $number_of_slices} {incr l} { 
  set counting_array($l) [expr $counting_array($l) * 1.0 / ( $total_number_frames * $slab_volume ) ]
  puts $fid "$l $counting_array($l)"
}
close $fid
mol delete top

exit
