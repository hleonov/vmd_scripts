# Secondary structure retension - Ugly script
# Computes percentage of a specific secondary structure.
# Ugly because it deals with the one-time ss calculation of VMD by writing and loading each frame. 

check_for_help "vmd -dispdev text -e ~/vmd_scripts/ss_ret_sel -args <pdb> <out> <trj-type> <skip> <trj> <H|C|T>"

set pdb			[lindex $argv 0];
set out			[lindex $argv 1];
set trj_type	[lindex $argv 2];
set skip    	[lindex $argv 3]; 
set trj 		[lindex $argv 4];
set ss			[lindex $argv 5];
set sel		[lreplace $argv 6]

mol load gro $pdb;
set fid [open $out w];

if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}

set prot [atomselect top "(resid 1 to 25) or (resid 39 to 63) or (resid 77 to 101) or (resid 115 to 139)" frame 0]
set ss_sel [atomselect top "((resid 1 to 25) or (resid 39 to 63) or (resid 77 to 101) or (resid 115 to 139)) and structure $ss" frame 0]

set number_of_frames [molinfo top get numframes]
for {set i 0} {$i < $number_of_frames} {incr i} {
	report $i 50
	#set prot [atomselect top "protein" frame $i]
	$prot frame $i
	$prot writepdb "tmp.pdb"
	mol load pdb "tmp.pdb"
	set ss_sel [atomselect top "((resid 1 to 25) or (resid 39 to 63) or (resid 77 to 101) or (resid 115 to 139)) and structure $ss" frame $i]
	mol delete top
	puts "$i [$ss_sel num] [$prot num]"
	set ret [expr 100*[$ss_sel num]/[$prot num]]
	puts $fid "$i\t$ret"
}

close $fid;
exit;
