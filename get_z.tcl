# Usage example:
# vmd -e get_z.tcl -args system_from_em.pdb system.xtc 10 resname AMA
source /Users/hleonov/vmd_scripts/utilities.tcl

check_for_help "vmdt -dispdev text -e ~/vmd_scripts/get_z.tcl <pdb> <xtc> <skip_frames> <selection: e.g resname AMA>";
set pdb		[lindex $argv 0];
set xtc 	[lindex $argv 1];
set skip_fr [lindex $argv 2]; 
set seltxt  [lreplace $argv 0 2];

set out "z_out.txt";
set fid [open $out w]
mol load pdb $pdb  
if {$skip_fr > 1} {
	animate read xtc $xtc skip $skip_fr waitfor all
} else {
	animate read xtc $xtc waitfor all
}

set number_of_frames [molinfo top get numframes]
set sel [atomselect top "$seltxt" frame 0];
for {set i 0} {$i < $number_of_frames} {incr i} {
	$sel frame $i
	set z [lindex [measure center $sel weight [$sel get mass]] 2]
	puts $fid "$z"
}

exit
