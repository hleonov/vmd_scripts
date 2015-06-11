# it runs on a single xtc and takes a while since it uses bigdcd to read one frame at a time
# the residue number are hard wired as is the selection to creat the cylinder around it

# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/relative_z.tcl -args after_6000.pdb system.xtc z_position_163.xls  z_position_164.xls 155 156
source ~arkini/data/scripts/Tcl/bigdcd.tcl
source ~arkini/data/scripts/Tcl/waitfor.tcl

# the command line arguments parsing
set pdb               [lindex $argv 0] ;
set xtc               [lindex $argv 1] ;
set file_a            [lindex $argv 2] ;
set file_b            [lindex $argv 3] ;
set residue_a         [lindex $argv 4] ;# 155
set residue_b         [lindex $argv 5] ;# 156

set slab_width 0.2
set total 350.0
set fid_a [open $file_a w]
set fid_b [open $file_b w]
# here I load the protein  
mol load pdb $pdb  
# run thru the trajectory
proc analyze {frame} {
  global slab_width total fid_a fid_b residue_a residue_b
  set all [atomselect top "all"]
  set ca [atomselect top "protein and type CA"]
  # first residue
  set sel [atomselect top "resid $residue_a and type OD1 OD2 and protein"]
  set z_system [lindex [measure center $all] 2]
  set z_sell [lindex [measure center $sel] 2]
  set diff_a [expr $z_sell - $z_system]
  set diff_a [expr $diff_a / $slab_width]
  set diff_a [expr ($total / 2) + $diff_a]
  puts $fid_a "$frame	 $diff_a"
  # second residue
  set sel [atomselect top "resid $residue_b and type OD1 OD2 and protein"]
  set z_system [lindex [measure center $all] 2]
  set z_sell [lindex [measure center $sel] 2]
  set diff_b [expr $z_sell - $z_system]
  set diff_b [expr $diff_b / $slab_width]
  set diff_b [expr ($total / 2) + $diff_b]
  puts $fid_b "$frame	 $diff_b"
  # printout every 50 frames
  if { $frame % 50 == 0} {
    puts "Analyzing frame: $frame $diff_a $diff_b"
  }
}
proc end {} {
  exit
}
waitfor $xtc end
bigdcd xtc analyze $xtc















