# (dx*dy / number of lipid in each leaflet)
source ~/vmd_scripts/utilities.tcl

#check_for_help "Usage: vmdt -e ~/vmd_scripts/area_per_lipid.tcl -args <gro> <trj-type <trj>";
set gro 		[lindex $argv 0]
set trj_type	[lindex $argv 1];
set trj 		[lindex $argv 2];

mol load gro $gro
animate read $trj_type $trj waitfor all
set sel [atomselect top "resname SOL"]
set fid [open "area_per_lipid.dat" w]
set lipids [atomselect top "type P"]
set leaflet [expr [$lipids num]/2]
set protein [atomselect top "protein"]
measure minmax $protein
for {set i 0} {$i <= [molinfo top get {numframes}]} {incr i} {

    $sel frame $i
	#set size [measure minmax $sel]
	set x_min [lindex [lindex [measure minmax $sel] 0] 0]
	set y_min [lindex [lindex [measure minmax $sel] 0] 1]
	set x_max [lindex [lindex [measure minmax $sel] 1] 0]
	set y_max [lindex [lindex [measure minmax $sel] 1] 1]
	
	set dx [expr $x_max - $x_min]
	set dy [expr $y_max - $y_min]
	set lipid_area [expr $dx * $dy / $leaflet]
	puts $fid "$i $lipid_area"
  
  }
  exit
#}
