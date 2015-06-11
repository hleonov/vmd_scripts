#This script moves amantadine to different starting positions in the channel. 
#starting from center of channel and moving outwards (ACE-n-terminal side)
#by dz Angstrom steps. or inward (C-terminal side) by -dz. both work until a certain min / max Z were acheived.
#Each position is built on top of the previous one (by deducting -dz from it). 
#The center of mass of AMN is being repositioned, not NAK. 
# run : vmd -dispdev text -e ./position_amn.tcl -args initial_system.pdb <dz> <z_trans =1|0> minZ maxZ drug-name


#USER NOTES: 1. which residues are saved if at all
#			 2. dz in command line is in absolute value
# 			 3. for the translation of amantadine only in Z - use z_trans=1
source /Users/hleonov/vmd_scripts/utilities.tcl
check_for_help "vmdt -dispdev text -e ./position_amn.tcl -args initial_system.pdb <dz> <z_trans =1|0> minZ maxZ drug-name"

mol load pdb [lindex $argv 0]
set in_dz 	 [lindex $argv 1]
set z_trans  [lindex $argv 2]
set min_Z 	 [lindex $argv 3]
set max_Z 	 [lindex $argv 4]
set drug_res [lindex $argv 5]

proc center_drug {} {
	global prot drug drug_res
	set v1 [measure center $prot]
	#set v2 [lindex [$amn_n get {x y z}] 0]
	set v2 [measure center $drug weight [$drug get mass]]
	#set v2 [measure center $drug]
	set vec_dist [vecsub $v1 $v2]; #vector from drug-->protein
	$drug moveby $vec_dist
	puts "Centering drug $drug_res in the protein"
}

# First position - center of channel
set prot [atomselect top "backbone"]
set drug [atomselect top "resname $drug_res"]
#set amn_n [atomselect top "type NAK"]

center_drug

#check 
#measure center $prot
#$amn_n get {x y z}

#save
set all [atomselect top "all"] 
#or 
#set all [atomselect top "all and not resid 199 170"]
set num 0
$all writepdb center_$num.pdb

#set low_d [expr $dz/2]
#set high_d [expr $dz+$low_d]

#dz is the new Z position here. not a delta value. 
proc move_drug {dz} {
	global z_trans drug num v1 v2 prot vec_dist
	
	if {$z_trans == 1} {
		$drug moveby [list 0 0 $dz]
	} else {
		set v2 [measure center $drug weight [$drug get mass]]
		#take next z position (v2_z + $dz) and average protein coordinates within epsilon=2A
		#to each side (on z axis)
		set low [expr [lindex $v2 2]+$dz-4]
		set high [expr [lindex $v2 2]+$dz+4]
		
		#if $dz is negative then $low and $high replace eachtother
		if {$low<$high} {
			set sel [atomselect top "(protein) and (z>=$low) and (z<=$high)"]
		} else {
			set sel [atomselect top "(protein) and (z>=$high) and (z<=$low)"]
		}
		set v1 [measure center $sel]
		set vec_dist [vecsub $v1 $v2]
		#lreplace <list> <start_i> <end_i> value (replace avg z by our dz)
		set vec_dist [lreplace $vec_dist 2 2 $dz] 
	
		$drug moveby $vec_dist
	}

	set all [atomselect top "all"] 
	#set all [atomselect top "all and not resid 189"] 
	 
	set pnum [expr $num+$dz]
	$all writepdb center_$pnum.pdb
}

# if we only need a short range, all negative, start with amn from max_Z (instead of 0) and decend to min_Z.
set init 0;
if {$max_Z < 0} {
	set init $max_Z;
	move_drug [expr $max_Z]
}
#progress by steps of dz on z-axis and adjust to pore center + save
if {$min_Z < 0} {
	for {set num $init} {$num>$min_Z} {set num [expr $num-$in_dz]} {
		move_drug [expr -1*$in_dz]
	}
}
center_drug

# if we only need a short range, all positive, start from min_Z and ascend to max_Z.
set init 0;
if {$min_Z > 0} {
	set init $min_Z;
	move_drug [expr $min_Z]
}
if {$max_Z>0} {
	for {set num $init} {$num<$max_Z} {set num [expr $num+$in_dz]} {
		move_drug $in_dz
	}
}
exit
