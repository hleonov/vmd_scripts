proc align {sel all mol} {
   for {set frame 0} {$frame < [molinfo $mol get numframes]} {incr frame} {
     $sel frame $frame
     $all frame $frame
     puts [measure center $sel]
     $all moveby [vecinvert [measure center $sel]]
   }
}

# Create selections
set all [atomselect top "all"]
set sel [atomselect top "protein and name CA"]
#
# Execute 
#
align $sel $ref $all top


