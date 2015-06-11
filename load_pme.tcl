
# set pdb [lindex $argv 0]
# set xtc [lindex $argv 1]
# set skip [lindex $argv 2]
# set res_sel [lindex $argv 3]; #"AMA HIS"
# set charges [lindex $argv 4]; #"charges_20ns.dx"
set pdb "V27A_with_charges.pdb"
set xtc "final_fit.xtc"
set skip 10
set res_sel "GLY"
set charges "new_charges_out.dx"
#[lreplace $argv 0 2]

#load system and color
mol load pdb $pdb
if {($skip>1)} {
	animate read xtc $xtc skip $skip waitfor all
} else {
	animate read xtc $xtc waitfor all
}

mol addrep top
mol modselect 0 top "protein"
mol modselect 1 top "resname $res_sel"
mol modstyle 0 top NewCartoon
mol modstyle 1 top Licorice
mol modcolor 0 top ColorID 2
mol modcolor 1 top Name


#set the charges from beta field and position at 0,0,0
set sel [atomselect top "all"]
$sel set charge [$sel get {beta}]
for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
	molinfo top set frame $i; 
	$sel moveby [vecinvert [measure center $sel]];
}

#activate PME - this gives weird color scheme
#set max [expr [molinfo top get numframes]-1]
#pmepot -sel $sel -frames 1:$max -dxfile charges_out.dx 

#if already calculated PME
mol new $charges type dx
