source ~/vmd_scripts/utilities.tcl

proc align {sel ref all mol} {
   for {set frame 0} {$frame < [molinfo $mol get numframes]} {incr frame} {
      $sel frame $frame
      $all frame $frame
      $all move [measure fit $sel $ref]
   }
}

#check_for_help "vmdt -dispdev text -e ~/vmd_scripts/load_for_pme.tcl <pdb> <xtc> <skip_frames>";

set pdb "V27A_with_charges.pdb";
set xtc "final_fit.xtc";
set skp 10;

mol load pdb $pdb
animate read xtc $xtc skip $skp waitfor all

set all [atomselect top "all"]
set sel [atomselect top "protein and type CA"]
set ref [atomselect top "protein and type CA" frame 0]

$all set charge [$all get {beta}]


for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
	molinfo top set frame $i; 
	$all moveby [vecinvert [measure center $all]]
	$all move [measure fit $sel $ref]
}



#align $sel $ref $all top

