# Secondary structure retension - not working!!
# Computes percentage of a specific secondary structure.

check_for_help "vmd -dispdev text -e ~/vmd_scripts/ss_retension -args <pdb> <out> <trj-type> <skip> <trj> <H|C|T>"
source "~/vmd_scripts/ss_utils.tcl"

set pdb			[lindex $argv 0];
set out			[lindex $argv 1];
set trj_type	[lindex $argv 2];
set skip    	[lindex $argv 3]; 
set trj 		[lindex $argv 4];
set ss			[lindex $argv 5];

mol load pdb $pdb;
set fid [open $out w];

if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}
start_sscache

set prot [atomselect top "protein" frame 0]
set ss_sel [atomselect top "protein and structure $ss" frame 0]

set number_of_frames [molinfo top get numframes]
for {set i 0} {$i < $number_of_frames} {incr i} {
	molinfo top set frame $i
	#set prot [atomselect top "protein" frame $i]
	$prot frame $i
	#vmd_calculate_structure top
	$ss_sel frame $i
#	$prot writepdb "tmp.pdb"
#	mol load pdb "tmp.pdb"
#	set ss_sel [atomselect top "protein and structure $ss" frame $i]
#	mol delete top

	set ret [expr 100*[$ss_sel num]/[$prot num]]
	puts $fid "$i\t$ret"
	puts "$i [$ss_sel num] [$prot num]"
	report $i 50
}
stop_sscache
close $fid;
exit;
