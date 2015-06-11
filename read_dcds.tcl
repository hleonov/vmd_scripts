# this script reads several dcds into the top molecule the form of 
# name$i.dcd, whereby $i runs from start to end, inclusive
# if $fit is equal to 1 it fits the Ca of all frames according to frame 0

proc read_dcds { name start end fit } {
  for {set i $start} {$i <= $end} {incr i} {
    set dcd_file "$name$i.dcd"
    animate read dcd $dcd_file waitfor all
  }
  if { $fit == 1 } {
    puts "Fitting all frames to inital frame"
    set all [atomselect top "all"]
    set sel [atomselect top "protein and name CA"]
    set ref [atomselect top "protein and name CA" frame 0]
    for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame} {  
       $sel frame $frame						      
       $all frame $frame						      
       $all move [measure fit $sel $ref]				      
    }									      
  }
}

