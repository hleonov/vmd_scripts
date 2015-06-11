# this script will operate on the top molecule
# 1. a pdb
# 2. a trajectory
# 3. a starting frame
# 4. an end frame
# 5. an index for atom j 
# 6. an index for atom j 
# 7. file name
#
# It will then generate an ascii list with the distance between the atoms (i,j) vs. frames.
# Usage:
# vmd -dispdev text -e get_z.tcl -args system_from_em.pdb system.xtc 0 2000 2389 67174 cg_k.xls


proc get_distance {pdb trajectory start end indexi indexj file} {
  set fid [open $file w]
  mol load pdb $pdb  
  animate read xtc $trajectory beg $start end $end waitfor all
  set number_of_frames [molinfo top get numframes]
  for {set i 0} {$i < $number_of_frames} {incr i} {
    # the index in vmd is zero based
    set j [expr $indexi - 1.0 ]
    set k [expr $indexj - 1.0 ]
    set sel_j [atomselect top "index $j" frame $i]
    set sel_k [atomselect top "index $k" frame $i]
    set distance  [veclength [vecsub [lindex [$sel_j get {x y z}] 0] [lindex [$sel_k get {x y z}] 0]]]
    puts $fid $distance
    puts $distance
  }
  animate delete beg 1
}
# run it as such
# vmd -dispdev text -e ~/data/vmd_scripts/get_distance.tcl -args last.pdb system.xtc 0 2000 2389 67174 k_z.xls
get_distance [lindex $argv 0] [lindex $argv 1] [lindex $argv 2] [lindex $argv 3] [lindex $argv 4] [lindex $argv 5] [lindex $argv 6] 
exit
