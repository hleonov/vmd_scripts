# This script computes the tilt angle of each helix from the bilayer normal. 
# Note that there are actually two bilayers normals - top or bottom leaflets.
# 

source ~hleonov/vmd_scripts/utilities.tcl

proc avg_tilt {frame resid_b resid_e n_top n_bot shift} {
	set sum_top 0
	set sum_bot 0
	for {set r $resid_b} {$r<=[expr $resid_e - $shift]} {incr r} {
		set first $r;
    	set last  [expr $r+$shift];
		set mytilt($r) [tilt $frame $first $last $n_top $n_bot]
		set sum_top [expr $sum_top + [lindex $mytilt($r) 0]]
		set sum_bot [expr $sum_bot + [lindex $mytilt($r) 1]]
    }
	
	return [list [expr $sum_top/[array size mytilt]] [expr $sum_bot/[array size mytilt]]]
}

proc tilt {frame first last n_top n_bot} {
  set start [atomselect top "resid $first to [expr $first + 3] and type CA" frame $frame]
  set end   [atomselect top "resid $last to[expr $last  - 3] and type CA" frame $frame]
  set start_coords [measure center $start]
  set end_coords   [measure center $end]
  set helix [vecsub $start_coords $end_coords]
  set angle_top [angle_between_two_vectors $helix $n_top]
  if { $angle_top > 90} {
    set angle_top [expr 180 - $angle_top]
  } 
  set angle_bot [angle_between_two_vectors $helix $n_bot]
  if { $angle_bot > 90} {
    set angle_bot [expr 180 - $angle_bot]
  } 
  #set angle_z [angle_between_two_vectors $helix {0 0 1}]
#  if { $angle_z > 90} {
#    set angle_z [expr 180 - $angle_z]
#  } 
  return [list $angle_top $angle_bot]
}

proc leaflet_normal {seltxt fr} {
	set I [pca $seltxt $fr]
	puts "I: $I"
	set sortI [lsort -command sort_by_z $I]
	set PCA [lreplace $sortI 2 2]
	set n [norm_to_plane [lindex $PCA 0] [lindex $PCA 1]]
	puts "normal: $n"
	return $n
}

proc def_helices {hx_bounds} {
	set helices [split $hx_bounds ":"];
	for {set j 0} {$j<[llength $helices]} {incr j} {
		set hx [lindex $helices $j]
		set hx_ends($j) [split $hx "-"];
	#	set sel [atomselect top "type CA and resid [lindex $hx_ends($j) 0] to [lindex $hx_ends($j) 1]" frame 0]
	#	$sel global
	#	set hx_sel($j) $sel
	#	puts "atoms in helix [$sel text] is [$hx_sel($j) num]"
	} 
	return [array get hx_ends]
}


check_for_help "vmdt -e ~/vmd_scripts/helix_tilt_trj.tcl -args <gro> <trj-type> <skip> <trj> <out> <helix-bound>"
set pdb			[lindex $argv 0];
set trj_type	[lindex $argv 1];
set skip    	[lindex $argv 2]; 
set trj 		[lindex $argv 3];
set out			[lindex $argv 4];
set hx_bounds 	[lreplace $argv 0 4]; 

#set out "tilt"
#set hx_bounds "1-25:26-50:51-75:76-100"


set shift 7

mol load pdb $pdb;

if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}

center_by_selection "protein and type CA" 0
set ref [atomselect top "protein and type CA" frame 0]

#select helices - possible improvement - no need for actual selection, only limits.
array set hx_sel [def_helices $hx_bounds]

#open output files - one for each helix
for {set j 0} {$j<[array size hx_sel]} {incr j} {
	set fid($j) [open "$out\_h$j\.dat" w];
	set avg($j) [open "avg_$out\_sh$shift\_h$j\.dat" w];
}

set number_of_frames [molinfo top get numframes]
for {set fr 0} {$fr < $number_of_frames} {incr fr} {

	report $fr 50
	
	#align to ref which was centered at 0 0 0
	set all [atomselect top "all" frame $fr]
	set ca [atomselect top "protein and type CA" frame $fr]
	$all move [measure fit $ca $ref]	
	
	#recompute membrane normal
	set n_top [leaflet_normal "type P and z>=0" $fr]
	set n_bot [leaflet_normal "type P and z<=0" $fr]
	
	#compute tilt
	for {set j 0} {$j<[array size hx_sel]} {incr j} {
		set tilt [tilt $fr [lindex $hx_sel($j) 0] [lindex $hx_sel($j) 1] $n_top $n_bot]
		set avg_tilt [avg_tilt $fr [lindex $hx_sel($j) 0] [lindex $hx_sel($j) 1] $n_top $n_bot $shift]
		puts $fid($j) "$fr $tilt"
		puts $avg($j) "$fr $avg_tilt"
	}
}

for {set j 0} {$j<[array size hx_sel]} {incr j} {
   close $fid($j)
   close $avg($j)
}

exit

