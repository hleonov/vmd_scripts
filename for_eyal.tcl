# the file name to write the reuslts in
set file "my_results.xls"
# the pdb
set pdb "protein.pdb"
# the trajectory in xtc format
set xtc "protein.xtc"
# the start point of the trajectory in ps
set start  0
# the end point of the trajectory in ps
# if you want to read all of it just chaneg the command below
set start  1000
# the index of the Ca atom
set i 100
# the index of atom O atom
set j 200

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
