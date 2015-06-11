# This will rotate the selection around a bond that is formed by atoms of index 1 and 2
proc rotate_bond { sel ind1 ind2 angle {molid top}} {
  set a1 [atomselect $molid "index $ind1"]
  set a2 [atomselect $molid "index $ind2"]
  set c1 [lindex [$a1 get {x y z}] 0]
  set c2 [lindex [$a2 get {x y z}] 0]
  # compute rotation matrix
  set rot_mat [trans bond $c1 $c2 $angle deg]
  # rotate selection
  $sel move $rot_mat
} 
