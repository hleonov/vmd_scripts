# This script moves a ligand to different starting positions in a protein (usually channel).
# starting from center of channel and moving outwards by dz Angstrom steps, or inward by (-dz), until a certain min / max Z were acheived.
# Each position is built on top of the previous one (by deducting/adding to it). 
# The center of mass of the ligand is being repositioned, not a specific atom.. 
# Run example: vmdt -e ~/vmd_scripts/position_drug.tcl -args m2+amn.pdb 1 1 -28 -20 AMA
# Usage : vmdt -e ./position_drug.tcl -args <initial.pdb> <dz> <z_trans =1|0> <minZ> <maxZ> <drug-name>
# Notes: initial PDB must contain the ligand (doesn't matter where).
#USER NOTES: 1. dz in command line is in absolute value
# 			 2. for the translation of ligand only in Z - use z_trans=1. 
#				Useful when ligand emerges into the solvent, and centering by protein coordinates is irrelevant.
#			 3. initial pdb must contain ligand, location doesn't matter.

source /Users/hleonov/vmd_scripts/utilities.tcl

proc center_drug {} {
	global prot drug drug_res
	set v1 [measure center $prot];							#v1 always represents COM of next location
	set v2 [measure center $drug weight [$drug get mass]];	#v2 always represents current drug COM
	set vec_dist [vecsub $v1 $v2]; 							#vec_dist is always the vector "drug-->protein"
	$drug moveby $vec_dist
	puts "Centering drug $drug_res in the protein"
}

proc move_drug {dz} {
	global z_trans drug num v1 v2 prot vec_dist
	
	if {$z_trans == 1} {
		$drug moveby [list 0 0 $dz]
	} else {
		set v2 [measure center $drug weight [$drug get mass]]
		#take next z position (v2_z + $dz) and average protein coordinates within epsilon=4A to each side (on z axis)
		set low [expr [lindex $v2 2]+$dz-4]
		set high [expr [lindex $v2 2]+$dz+4]
		
		#select protein atoms around the next z position - so that the drug is positioned in their center.
		#if $dz is negative then $low and $high replace eachother
		if {$low<$high} {
			set sel [atomselect top "(protein) and (z>=$low) and (z<=$high)"]
		} else {
			set sel [atomselect top "(protein) and (z>=$high) and (z<=$low)"]
		}
		set v1 [measure center $sel]
		set vec_dist [vecsub $v1 $v2]
		set vec_dist [lreplace $vec_dist 2 2 $dz]; 			#correct to moveby dz precisely 
	
		$drug moveby $vec_dist
	}

	set all [atomselect top "all"] 
	#set all [atomselect top "all and not resid 189"] 
	 
	set pnum [expr $num+$dz]
	$all writepdb center_$pnum.pdb
}

check_for_help "vmdt -dispdev text -e ./position_drug.tcl -args <initial.pdb> <dz> <z_trans =1|0> <minZ> <maxZ> <drug-name>"

mol load pdb [lindex $argv 0]
set in_dz 	 [lindex $argv 1]
set z_trans  [lindex $argv 2]
set min_Z 	 [lindex $argv 3]
set max_Z 	 [lindex $argv 4]
set drug_res [lindex $argv 5]

# First position - center of channel (center_0).
set prot [atomselect top "backbone"]
set drug [atomselect top "resname $drug_res"]

center_drug

#save center 0
set all [atomselect top "all"] 
set num 0
$all writepdb center_$num.pdb

# if we only need a short range, all negative, start with ligand at max_Z (instead of 0)
set init 0;
if {$max_Z < 0} {
	set init $max_Z;
	move_drug [expr $max_Z]
}

#decend by steps of dz on z-axis and adjust to pore center + save, until min_Z is reached.
if {$min_Z < 0} {
	for {set num $init} {$num>$min_Z} {set num [expr $num-$in_dz]} {
		move_drug [expr -1*$in_dz]
	}
}
center_drug

# if we only need a short range, all positive, start from min_Z.
set init 0;
if {$min_Z > 0} {
	set init $min_Z;
	move_drug [expr $min_Z]
}
#ascend by steps of dz to max_Z
if {$max_Z>0} {
	for {set num $init} {$num<$max_Z} {set num [expr $num+$in_dz]} {
		move_drug $in_dz
	}
}
exit
