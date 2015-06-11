# (dx*dy / number of lipid in each leaflet)
#source ~/vmd_scripts/shy_utilities.tcl

#check_for_help "Usage: vmdt -e ~/vmd_scripts/area_per_lipid.tcl -args <gro> <xtc>[-h]";
set gro [lindex $argv 0]
set trjtype [lindex $argv 1]
set trj 	[lindex $argv 2]
set mc_n	[lindex $argv 3]

mol load gro $gro
animate read $trjtype $trj waitfor all
#  set sel [atomselect top "resname SOL"]
set fid [open "area_per_lipid.dat" w]
for {set i 0} {$i <= [molinfo top get {numframes}]} {incr i} {
	#    $sel frame $i
	
	set lipid_area [expr $dx * $dy / 60.0]
	puts $fid "$i $lipid_area"
  
  }
  exit
#}
