# this will rotate and tilt the helix as requested
source /Users/stan/Bin/utilities.tcl

if {[lindex $argv 0] == "-rf"} {
  set pdb       [lindex $argv 1];
  set out       [lindex $argv 2];
  set by_tilt   [lindex $argv 3];
  set by_rotate [lindex $argv 4];
} else {
  set pdb       [lindex $argv 0];
  set out       [lindex $argv 1];
  set by_tilt   [lindex $argv 2];
  set by_rotate [lindex $argv 3];
}



mol load pdb $pdb

set all [atomselect top "all"];
set mat [rotation_x $by_tilt];
$all move $mat;
set mat2 [rotation_helix $by_rotate 4 0];
$all move $mat2;
$all writepdb [out_file $out];

exit











