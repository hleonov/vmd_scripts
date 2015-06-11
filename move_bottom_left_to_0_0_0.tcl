# will move the system such that th ebottom left corner is at 0 0 0 
# set box size
# usage:
# vmdt -e  /Users/hleonov/vmd_scripts/move_bottom_left_to_0_0_0.tcl -args system.pdb system_moved.pdb

source /Users/hleonov/vmd_scripts/shy_utilities.tcl
set in  [lindex $argv 0];
set out [lindex $argv 1];
set out [out_file $out]

mol load pdb $in
set sel [atomselect top "all"]
$sel moveby [vecinv [lindex [measure minmax $sel] 0]]
set water [atomselect top "resname SOL HOH"]
set sel [atomselect top "all"]
set size [measure minmax $sel]
molinfo top set {a b c} [lindex $size 1]
$sel writepdb $out
exit
