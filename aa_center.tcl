source ~/vmd_scripts/utilities.tcl

# Reads a gro file and following trajectories, and computes the Z distribution of each residue along the pore.
# The average is given on all 4 helices together. 

check_for_help "vmd -dispdev text -e ~/vmd_scripts/aa_center.tcl -args <minZ> <maxZ> <Dz> <incl. bb 1/0> <file_pre> <file_suf>"

set minZ		[lindex $argv 0]
set maxZ		[lindex $argv 1]
set dz			[lindex $argv 2]	
set bb			[lindex $argv 3]
set file_pre	[lindex $argv 4]
set file_suf	[lindex $argv 5]

set resi_begin 1
set resi_end   25

for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set z_center_arr($resi) 0
}
set fid [open aa_center.dat w]

puts "before for minZ $minZ maxZ $maxZ"
for {set coord $minZ} {$coord <= $maxZ} {}  {
	puts "current: $coord"
    set coord_num $coord
	set filename "$file_pre\_$coord_num$file_suf.pdb"
	puts "$filename"
	mol load pdb $filename
  	for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set id2 [expr $resi + $resi_end]
	set id3 [expr $id2 + $resi_end]
	set id4 [expr $id3 + $resi_end]
	if {$bb == 1} {
		set all_amino_sel [atomselect top "resid $resi $id2 $id3 $id4" ]
	} else {
		set all_amino_sel [atomselect top "resid $resi $id2 $id3 $id4 and not backbone" ]
	}
	set z_center_arr($resi) [expr $z_center_arr($resi) + [lindex [measure center $all_amino_sel] 2]]
#    mol delete top
  }
  set coord [expr $coord + $dz]
  	 
}

for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set aa [atomselect top "type CA and resid $resi" ]
	set name [$aa get resname]
	set z_center [expr $z_center_arr($resi) / [expr ($maxZ - $minZ + 1)/$dz]]
	puts $fid "$resi $name $z_center"
}
exit

#puts "$resi $id2 $id3 $id4\n"
	#foreach pres [$all_amino_sel get {name resid resname chain z}] {
	#	puts "$pres"
	#}
