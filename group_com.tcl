source ~/vmd_scripts/utilities.tcl

# Reads a pdb file and following trajectory, plus a selection text
# computes and writes the Center Of Mass for selection for each loaded frame.
# Usage example:
# vmd -dispdev text -e ~/vmd_scripts/group_com.tcl -args center_-3.pdb 1 from_umb_center_-3.xtc out.dat resname AMA

check_for_help "vmd -dispdev text -e ~/vmd_scripts/group_com.tcl -args <pdb> <skip> <trj> <out> <seltext>"

set pdb			[lindex $argv 0];
set skip    	[lindex $argv 1]; 
set trj 		[lindex $argv 2];
set out			[lindex $argv 3];
set seltext		[lreplace $argv 0 3];

mol load pdb $pdb;
set group [atomselect top "$seltext"];
set fid [open $out w]

#read the trajectory file
if {($skip>1)} {
	animate read xtc $trj skip $skip waitfor all
} else {
	animate read xtc $trj waitfor all
}

for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
      $group frame $i
      set com [measure center $group]
	  set com_z [lindex $com 2]
      puts $fid "$i $com_z"
	  flush $fid; # not really needed
}

close $fid;
unset group
exit


