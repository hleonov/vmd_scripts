# Input: PDB filename, resname, {x y z} shift, outname 
# Output: new PDB with resname moved by xyz,

source /Users/hleonov/vmd_scripts/utilities.tcl
check_for_help "vmdt -dispdev text -e ./relocate_drug.tcl -args <initial.pdb> <resname> <x y z> <out>"
mol load pdb [lindex $argv 0]
set move_res [lindex $argv 1]
set out 	 [lindex $argv 2]
set x		 [lindex $argv 3]
set y		 [lindex $argv 4]
set z		 [lindex $argv 5]

set to_move [atomselect top "resname $move_res"]
set c [measure center $to_move]
$to_move moveby [list $x $y $z]
#update $to_move
set c [measure center $to_move]

set all [atomselect top "all"]
$all writepdb $out
exit
