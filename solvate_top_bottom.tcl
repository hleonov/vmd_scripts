# this adds a layer of water to the top and bottom of the membrane
# the inputs are the thickness of the top  and bottom layers
# the files are called top_layer.pdb and bottom_layer.pdb
package require solvate
proc solvate_top_bottom {top bottom} {
  set sel [atomselect top "all"]
  $sel moveby [vecinvert [measure center $sel]]
  set a [expr [molinfo top get {a}] / 2]
  set b [expr [molinfo top get {b}] / 2]
  set c [expr [molinfo top get {c}] / 2]
  # the top
  if {$top>0} {
	  solvate -o top_layer -minmax [list [list [expr 0 - $a] [expr 0 - $b] [expr 0 + $c]] [list [expr 0 + $a] [expr 0 + $b] [expr 0 + $c + $top]]]
  }
  # the bottom
  if {$bottom>0} {
	  solvate -o bottom_layer -minmax [list [list [expr 0 - $a] [expr 0 - $b] [expr 0 - $c - $bottom]] [list [expr 0 + $a] [expr 0 + $b] [expr 0 - $c]]]
  }
  # these have it so that it isnt; centered so you don;t need to resave your original pdb
  # the top
  #solvate -o top_layer -minmax [list [list [expr 0 - $a + $a] [expr 0 - $b + $b] [expr 0 + $c + $c]] [list [expr 0 + $a + $a] [expr 0 + $b + $b] [expr 0 + $c + $top + $c]]]
  # the bottom
  #solvate -o bottom_layer -minmax [list [list [expr 0 - $a + $a] [expr 0 - $b + $b] [expr 0 - $c - $bottom + $c]] [list [expr 0 + $a + $a] [expr 0 + $b + $b] [expr 0 - $c + $c]]]
}




