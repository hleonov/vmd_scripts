# this script will operate on the top molecule
# 1. a pdb
# 2. a trajectory
# 3. a starting frame
# 4. an end frame
# 5. an index for the atom
# 6. file name
#
# It will then generate an ascii list with the z coordinate of the atom.
# is will fix all of the frames acording to the innitial frame
# Usage:
# vmd -dispdev text -e get_z.tcl -args system_from_em.pdb system.xtc 0 2000 56491 na_z.xls


proc get_z {pdb trajectory start end index file} {
  set fid [open $file w]
  mol load pdb $pdb  
  animate read xtc $trajectory beg $start end $end waitfor all
  set number_of_frames [molinfo top get numframes]
  set ref [atomselect top "protein and name CA" frame 0]
  for {set i 0} {$i < $number_of_frames} {incr i} {
    # Align and move
    set all [atomselect top "all" frame $i]
    set ca [atomselect top "protein and name CA" frame $i]
    $all move [measure fit $ca $ref]
    # the index in vmd is zero based
    set n [expr $index - 1.0 ]
    set sel [atomselect top "index $n" frame $i]
    puts $fid [$sel get {z}]
    # puts [$sel get {z}]
  }
  animate delete beg 1
}
#get_z after_12000.pdb system.xtc 0 100 67174 k_z.xls
# run it as such
# vmd -dispdev text -e ~/data/vmd_scripts/get_z.tcl -args after_12000.pdb system.xtc 0 500 67174 k_z.xls
get_z [lindex $argv 0] [lindex $argv 1] [lindex $argv 2] [lindex $argv 3] [lindex $argv 4] [lindex $argv 5] 
exit
