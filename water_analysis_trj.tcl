# This script computes a histogram of water throughout the channel. It loads each frame at a time. 
# In contrast to my script, it doesn't use a fit between the protein trajectories, but moves them all to [0 0 $range]
# Range is a given input instead of computed. 
# 

source /Users/hleonov/vmd_scripts/Tcl/bigdcd.tcl
source /Users/hleonov/vmd_scripts/Tcl/waitfor.tcl
source /Users/hleonov/vmd_scripts/utilities.tcl

check_for_help "Usage: vmd -dispdev text -e /Users/hleonov/vmd_scripts/water_analysis_trj.tcl -args system_from_PR0_com.gro after_fit.trr trr results.xls 50 35.0 3"

# the command line arguments parsing
set pdb               [lindex $argv 0] ;
set trj               [lindex $argv 1] ;
set trj_type          [lindex $argv 2] ;
set file			  [lindex $argv 3] ;	
set number_of_slices  [lindex $argv 4] ;
set range             [lindex $argv 5] ;
set radius            [lindex $argv 6] ;

# here I load the protein  
mol load pdb $pdb 
# the ouput file
set out [open [out_file $file] w];# check if the output file exists, add -rf in the command line to supress it
# define some basic parameters
set total_number_frames 0
# some simple calculations
set slab_width [expr 2 * $range * 1.0 / $number_of_slices ]
set r2 [expr $radius*$radius];
# initialze the counting array
for {set i 0} {$i <= $number_of_slices + 1} {incr i} {
  set counting_array($i) 0
}
# run thru the trajectory
proc analyze {frame} {
  global counting_array r2 range slab_width total_number_frames
  # the following centers the system around the selection and then shifts is by the vertical range
  # this is just for ease of computation. # note that the cylinder is now centered at the selection at EACH frame.
  # for the 1st residue
  set all [atomselect top "all"]
  set sel [atomselect top "protein and type CA"]
  $all moveby [vecinvert [measure center $sel]]
  $all moveby [list 0 0 $range]
  set water_sel [atomselect top "name OW and x^2 + y^2 <= $r2 and (z >= 0 ) and (z <= 2* $range )"]
  # now count the number of atoms in each histogram bin
  foreach z [$water_sel get {z}] {     
    set slice [ expr floor( [expr $z / $slab_width ] ) + 1 ] 
    set slice [ expr int($slice)] 
    incr counting_array($slice)
  }
  # this makes a nice printout every 50 frames
  if {[expr $total_number_frames / 50.0] == int([expr $total_number_frames / 50.0])} {
   puts "Analyzing frame: $total_number_frames"
  }
  incr total_number_frames
}

proc end {} {
  global total_number_frames r2 slab_width out number_of_slices counting_array
  set slab_volume [expr $r2 * 3.14159265358979 * $slab_width]
  for {set l 1} {$l < $number_of_slices} {incr l} { 
    set counting_array($l) [expr $counting_array($l) * 1.0 / ( $total_number_frames * $slab_volume ) ]
    puts $out "$l\t $counting_array($l)"
  }
  close $out
  exit
}

waitfor $trj end
bigdcd $trj_type analyze $trj

