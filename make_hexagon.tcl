# this sciprt will take a mambrane system and make a hexgon of it.
# it assumes that the membrane plane lies in the xy plane

# first center the molecule at 0 0 0 according to the water.
set water [atomselect top "water"]
set all [atomselect top "all"]
$all moveby [vecinvert [measure center $water]]
# now grab the x and y dimansion of the box again according to the water
set x [expr ([lindex [lindex [measure minmax $water] 1] 0] - [lindex [lindex [measure minmax $water] 0] 0]) / 2]
set y [expr ([lindex [lindex [measure minmax $water] 1] 1] - [lindex [lindex [measure minmax $water] 0] 1]) / 2]
# I now assuem the box is rectengular so lets average x and y to a singel parameter w
set w [expr ($x + $y) / 2]
# predefine the square root of 3
set s [expr sqrt ( 3 )]	
# select the top right triange
set top_right    "( y + ($s * x) - ($w * $s) > 0)"
set top_left     "( y - ($s * x) - ($w * $s) > 0)"
set bottom_right "(-y + ($s * x) - ($w * $s) > 0)"
set bottom_left  "(-y - ($s * x) - ($w * $s) > 0)"
set sel [atomselect top "same residue as (not ($top_right or $top_left or $bottom_right or $bottom_left))"]
$sel writepdb hexagon.pdb
