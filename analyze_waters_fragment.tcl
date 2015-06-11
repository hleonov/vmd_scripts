# this script should run through a trajectory and generate a histogram with waters
# the outcome is waters per A cubed. Bulck water density is 0.033.

# it runs on a single xtc and takes a while since it uses bigdcd to read one frame at a time
# the residue number are hard wired as is the selection to creat the cylinder around it
# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/analyze_waters_fragment.tcl -args after_6000.pdb system.xtc results.xls 50 35 7 500

source /Users/hleonov/vmd_scripts/Tcl/bigdcd.tcl
source /Users/hleonov/vmd_scripts/Tcl/waitfor.tcl
source /Users/hleonov/vmd_scripts/utilities.tcl

check_for_help "Usage: vmd -dispdev text -e ~/vmd_scripts/analyze_waters_fragment.tcl -args after_6000.gro system.xtc xtc results 50 35 7 500\n
				or\n
				vmdt -e ~/vmd_scripts/analyze_waters_fragment.tcl -args after_6000.gro system.xtc xtc results 50 35 7 500"

# the command line arguments parsing
set pdb               [lindex $argv 0] ;# system_from_em_free.pdb
set trj               [lindex $argv 1] ;# system.xtc
set trj_type		  [lindex $argv 2] ;# trr or xtc
set out	              [lindex $argv 3] ;# results
set number_of_slices  [lindex $argv 4] ;# 50
set range             [lindex $argv 5] ;# 35
set radius            [lindex $argv 6] ;# 9
set frame_stride      [lindex $argv 7] ;# 500

# here I load the protein  
mol load pdb $pdb  
# define some basic parameters
set total_number_frames 0
set stride_counter 0
# some simple calculations
set slab_width [expr 2 * $range * 1.0 / $number_of_slices ]
set r2 [expr $radius*$radius]
# initialze the counting array
for {set i 0} {$i <= $number_of_slices + 1} {incr i} {
  set counting_array($i) 0
}

proc write_out {} {
  global counting_array r2 range slab_width total_number_frames frame_stride out number_of_slices
  # open the file handle
  #set fid_b [open "${file_b}_${total_number_frames}.dat" w]
  set fid_out [open "${out}_${total_number_frames}.xls" w]
  # do the calcualtions are write-out
  set slab_volume [expr $r2 * 3.14159265358979 * $slab_width]
  for {set l 1} {$l < $number_of_slices} {incr l} { 
    set counting_array($l) [expr $counting_array($l) * 1.0 / ( $frame_stride * $slab_volume ) ]
    puts $fid_out "$l $counting_array($l)"
  }
  
  # close the file handles
  close $fid_out
  
  # zero the parameters
  for {set i 0} {$i <= $number_of_slices + 1} {incr i} {
    set counting_array($i) 0
  }

}
# run thru the trajectory
proc analyze {frame} {
  global counting_array r2 range slab_width total_number_frames frame_stride
  # the following centers the system around the selection and then shifts is by the vertical range
  # this is just for ease of computation. # note that the cylinder is now centered at the selection at EACH frame.
  # for the 1st residue
  set all [atomselect top "all"]
  set sel [atomselect top "protein and type CA"]
  $all moveby [vecinvert [measure center $sel]]
  # now move so that the box is centered such that its z is at 0
  $all moveby [vecinvert [list 0 0 [lindex [measure center $all] 2 ]]]
  # now move so that the box is centered such that its z is at - $range
  $all moveby [list 0 0 $range]
  # now the real selection
  set water_sel [atomselect top "name OW and (x^2 + y^2) <= $r2 and (z >= 0 ) and (z <= 2* $range )"]
#  set water_sel [atomselect top "name OW and (x^2 + y^2) <= $r2"]
#  puts "[$water_sel num]"
  # now count the number of atoms in each histogram bin
  foreach z [$water_sel get {z}] {     
    set slice [ expr floor( [expr $z / $slab_width ] ) + 1 ] 
    set slice [ expr int($slice)] 
    incr counting_array($slice)
  }
  incr total_number_frames
  # printout every 50 frames
  if { $total_number_frames % 50 == 0} {
    puts "Analyzing frame: $total_number_frames"
  }
  # this indicates that a batch of frames is ready to be processed
  if { $total_number_frames % $frame_stride == 0} {
    write_out
  }
}
proc end {} {
  
  exit
}
waitfor $trj end
bigdcd $trj_type analyze $trj
