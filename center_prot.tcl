set number_of_frames [molinfo top get numframes]
set start [atomselect top "all" frame 0]
set protein [atomselect top "protein or resname ACE NAC" frame 0]
set vecinv [vecinvert [measure center $protein]]
$start moveby [vecinvert [measure center $protein]]

set ref [atomselect top "protein and name CA" frame 0]
for {set k 1} {$k < $number_of_frames} {incr k} { 
  set all [atomselect top "all" frame $k]
  set ca [atomselect top "protein and name CA" frame $k]
  $all move [measure fit $ca $ref]
}
