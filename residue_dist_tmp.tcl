# This script is run like this: 
#   vmd -dispdev text -e ~/vmd_scripts/residue_dist.tcl -args system_from_PR0_com.gro trr 10 final_fit.trr his_dist.txt 16 41 66 9
# Where the last numbers represent resid numbers for residues for which to measure all-to-all distances.
# The output is a text file in which the first line is the legend, then the first column is frames, and the rest are distances


source ~/vmd_scripts/utilities.tcl

proc make_selections {index_list gsel} {
	global for_min, a1_sel;
	upvar $gsel sel
	set j 0;
	#select each residue in seperate
	foreach rsel $index_list {
		set j [expr $j+1];
		set sel($j) [atomselect top "resid $rsel and not backbone and noh" frame 0];
		#set id($j) "$rsel";

		#for min dist 
		set for_min($j) [atomselect top "resid $rsel" frame 0];
#		set sel1 [$r1 list];
#		set k 0;
#		foreach atom1 $sel1 {
#			set a1_sel($j,$k) [atomselect top "index $atom1" frame 0];
#			incr k;
#		}
	}
}
proc get_residue_dist {index_list fid mid gsel fr} {
	global for_min, a1_sel, sel;
	puts " the number is:  [$sel(1) num]"
	
	
	
	set j 0;
	#select each residue in seperate and calculate COM
	foreach rsel $index_list {
		set j [expr $j+1];
		
		$sel($j) frame $fr;
		puts " the number is:  [num $sel($j)]"
	#	set sel [atomselect top "resid $rsel and not backbone and noh" frame $fr];
		set v($j) [measure center $sel($j)];
		set id($j) "$rsel";
	}
	#measure distance between each COM
	puts -nonewline $fid "$fr\t";
	puts -nonewline $mid "$fr\t";
	for {set k1 1} {$k1 < $j} {incr k1} {
		for {set k2 $k1} {$k2 <= $j} {incr k2} {
			if {$k1 != $k2} {
#				set mind($k1,$k2) [get_min_residue_dist $id($k1) $id($k2) $fr];
#				set mind($k1,$k2) [get_min_residue_dist $k1 $k2 $fr $j];
				set dist($k1,$k2) [veclength [vecsub $v($k1) $v($k2) ]];
				puts -nonewline $fid "[format %8.3f $dist($k1,$k2)]\t";
				puts -nonewline $mid "[format %8.3f $mind($k1,$k2)]\t";
			}
		}
	}
	puts $fid "";
	puts $mid "";
}

proc get_min_residue_dist {r1 r2 fr j} {
	global for_min, a1_sel, sel;
	set sel1 [[atomselect top "resid $r1" frame $fr] list];
	set sel2 [[atomselect top "resid $r2" frame $fr] list];
	set min_dist 1000;
	
	#$for_min($k1) frame $fr;
#	$for_min($k2) frame $fr;
#	set sel1 [$r1 list];
#	set sel2 [$r2 list];
	#set b1 0;
	#set b2 0;
	foreach atom1 $sel1 {
		foreach atom2 $sel2 {
#			$a1_sel($j,$b1) frame $fr;
#		    $a2_sel($j,$b2) frame $fr;
#			set coord1 [lindex [$a1_sel($j,$b1) get {x y z}] 0];
#			set coord2 [lindex [$a1_sel($j,$b2) get {x y z}] 0];
			set a1_sel [atomselect top "index $atom1" frame $fr];
			set a2_sel [atomselect top "index $atom2" frame $fr];
			set coord1 [lindex [$a1_sel get {x y z}] 0];
			set coord2 [lindex [$a2_sel get {x y z}] 0];
			set a1a2_dist [veclength [vecsub $coord1 $coord2]];
			if {$a1a2_dist < $min_dist} {
				set min_dist  $a1a2_dist;
			}
		#	incr b2;
		}
		#incr b1;
	}
	return $min_dist;
}

proc get_pair_residue_dist {r1 r2 fr} {
	set sel1 [atomselect top "resid $r1 and not backbone and noh" frame $fr];
	set sel2 [atomselect top "resid $r2 and not backbone and noh" frame $fr];
	set c1 [measure center $sel1];
	set c2 [measure center $sel2];
	return [veclength [vecsub $c1 $c2 ]]
}

proc print_names {index_list fid mid} {
	# Make selections to get names
	set j 0;
	foreach rsel $index_list {
		set j [expr $j+1];
		set sel [atomselect top "resid $rsel"];
		set rname [lsort -unique [$sel get resname]]
		set name($j) [concat $rname$rsel];	
	}
	puts -nonewline $fid "# frame\t";
	puts -nonewline $mid "# frame\t";
	# print pairs of (resname_resid)<-->(resname_resid)
	for {set k1 1} {$k1 < $j} {incr k1} {
		for {set k2 $k1} {$k2 <= $j} {incr k2} {
			if {$k1 != $k2} {
				puts -nonewline $fid "$name($k1)\_$name($k2)\t"
				puts -nonewline $mid "$name($k1)\_$name($k2)\t"
			}
		}
	}
	puts $fid "";
	puts $mid "";
}
check_for_help "vmd -dispdev text -e ~/vmd_scripts/residue_dist.tcl -args <pdb> <trj-type> <skip> <trj> <out1> <out_min> <resid-list>"

set pdb			[lindex $argv 0];
set trj_type	[lindex $argv 1];
set skip    	[lindex $argv 2]; 
set trj 		[lindex $argv 3];
set out			[lindex $argv 4];
set out_min		[lindex $argv 5];
set res_sel 	[lreplace $argv 0 5]; 

for {set i 1} {$i <= 4} {incr i} {
  set gsel($i) 0;
}

mol load gro $pdb;
set fid [open [out_file $out] w];
set mid [open [out_file $out_min] w];
if {($skip>1)} {
	#animate read $trj_type $trj skip $skip waitfor all
	animate read $trj_type $trj beg 0 end 100 waitfor all
} else {
	animate read $trj_type $trj waitfor all
}
print_names $res_sel $fid $mid;
make_selections $res_sel $gsel;
set number_of_frames [molinfo top get numframes]
for {set i 0} {$i < $number_of_frames} {incr i} {
	report $i 50
	get_residue_dist $res_sel $fid $mid $gsel $i;
}

close $fid;
close $mid;
exit;
