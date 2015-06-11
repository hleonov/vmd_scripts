proc get_min_residue_dist {r1 r2 fr} {
	set sel1 [[atomselect top "resid $r1" frame $fr] list];
	set sel2 [[atomselect top "resid $r2" frame $fr] list];
	set min_dist 1000;
	foreach atom1 $sel1 {
		foreach atom2 $sel2 {
			set a1_sel [atomselect top "index $atom1" frame $fr];
			set a2_sel [atomselect top "index $atom2" frame $fr];
			set coord1 [lindex [$a1_sel get {x y z}] 0];
			set coord2 [lindex [$a2_sel get {x y z}] 0];
			set a1a2_dist [veclength [vecsub $coord1 $coord2]];
			if {$a1a2_dist < $min_dist} {
				set min_dist  $a1a2_dist;
			}
		}
	}
	return $min_dist;
}

proc read_list {list_file} {
	global pair;
	set i 0;
	set file_id [open $list_file r];
	#gets with 2 arguments reads into $line, and returns length, -1 for EOF
	while {[gets $file_id line] >= 0} {
		#split line into a list that can be accessed with lindex
		set pair($i) [split $line " "]; 
		incr i;		
	}
}


set pdb			[lindex $argv 0];
set plist 		[lindex $argv 1];
set out			[lindex $argv 2];
set trj_type	[lindex $argv 3];
set skip    	[lindex $argv 4]; 
set trj 		[lindex $argv 5];

mol load gro $pdb;
read_list $plist;
set fid [open $out w];

if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}

set number_of_frames [molinfo top get numframes]
for {set i 0} {$i < $number_of_frames} {incr i} {
	report $i 50
	set size [array size pair]
	for {set j 0} {$j<$size} {incr j} {
		get_min_residue_dist [lindex $pair($j) 0] [lindex $pair($j) 1] $i;
	}
}

close $fid;
exit;
