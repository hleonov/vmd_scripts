# This script is run like this: 
#   vmd -dispdev text -e ~/vmd_scripts/residue_dist.tcl -args system_from_PR0_com.gro trr 10 final_fit.trr his_dist.txt 16 41 66 9
# list dihedrals to compute: 4 atoms per line

source ~/vmd_scripts/utilities.tcl

proc compute_dihedrals {fr} {
	global atoms fid;
	set size [array size atoms];
	puts -nonewline $fid "$fr\t"
	for {set i 0} {$i < $size} {incr i} {
		set angle [dihedral [lindex $atoms($i) 0] [lindex $atoms($i) 1] [lindex $atoms($i) 2] [lindex $atoms($i) 3] $fr]
		puts -nonewline $fid "[format %8.3f $angle]\t";
	}
	puts $fid "";

}
proc read_list {list_file} {
	global atoms;
	set i 0;
	set file_id [open $list_file r];
	#gets with 2 arguments reads into $line, and returns length, -1 for EOF
	while {[gets $file_id line] >= 0} {
		#split line into a list that can be accessed with lindex
		set atoms($i) [split $line " "]; 
		incr i;		
	}

}

check_for_help "vmd -dispdev text -e ~/vmd_scripts/compute_dihedrals -args <pdb> <list> <out> <trj-type> <skip> <trj>"

set pdb			[lindex $argv 0];
set list_a 		[lindex $argv 1];
set out			[lindex $argv 2];
set trj_type	[lindex $argv 3];
set skip    	[lindex $argv 4]; 
set trj 		[lindex $argv 5];

mol load pdb $pdb;
read_list $list_a;
set fid [open $out w];

if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}

set number_of_frames [molinfo top get numframes]
for {set i 0} {$i < $number_of_frames} {incr i} {
	report $i 50
	compute_dihedrals $i;
}

close $fid;
exit;
