# calculate the condutnace events in the protein
# This works by defining the system into 3 segments:
# a. above the channel (z>0)
# b. in the channel    
# c. below the channel (z<0)
# any water can have 5 states from -2 to +2, as follows:
# -2 is below the channel
# -1 is in the channel, but before was below (i.e. it came in from below)
#  0 in the channel
#  1 is in the channel, but before was above it (i.e. it came in from above)
#  2 is above the channel
#
# In the beinging a water molecule above the channel is designated as "2"
#                 a water molecule below the channel is designated as "-2"
#                 a water molecule in    the channel is designated as "0"
# Keep 3 lists: the current state of water "state"
#		the 2nd list is just a list of all atoms with their indecies: index
#		the 3rd list is just a list of all atoms with their indecies: resid
# at every state form the new according to the position (i.e. it should only have values of  above, in, below
# then update the new list according to the old following these criteria:

# above
# before   should be
#   2	      2        	
#   1	      2
#   0	      2
#  -1	      EVENT UP 2
#  -2	      2

# in
# before   should be
#   2	      1
#   1	      1
#   0	      0
#  -1	     -1
#  -2	     -1

# below
# before   should be
#   2	     -2
#   1	     EVENT DOWN -2
#   0	     -2
#  -1	     -2
#  -2	     -2

source /Users/hleonov/vmd_scripts/Tcl/bigdcd.tcl
source /Users/hleonov/vmd_scripts/Tcl/waitfor.tcl
source /Users/hleonov/vmd_scripts/utilities.tcl

# the -rf needs to be the first argument due to the ability to grab many xtcs
check_for_help "Usage: vmd -dispdev text -e ~/vmd_scripts/conductance.tcl -args init_trj.pdb cond.log conductance.xls final_fit.xtc xtc filter"

set pdb  [lindex $argv 0];
set logf  [lindex $argv 1];
set outf  [lindex $argv 2];
set trj  [lindex $argv 3];
set trj_type [lindex $argv 4];
#set filter	[lreplace $argv 0 4]
set filter "resname TRP VAL"
set dt 10; # the timestep of each frame in ps
set radius 8; # the radius from the selectivity filter center

set dtop 0;
set dbot 0;
mol load pdb $pdb  
set out [open $outf w];
set log [open $logf w];
# get the number of waters
#set k [atomselect top "index 9665"]
set k [atomselect top "name OW and resname SOL"]
#set k [atomselect top "resname K"]
set n [$k num]
# the atom indices of the ions
set index [$k get index]
set resid [$k get resid]
# the selectivity filter 
#resid $filter
set sf [atomselect top "$filter"] 
# its edges
set top    [lindex [lindex [measure minmax $sf] 1] 2]
set bottom [lindex [lindex [measure minmax $sf] 0] 2]
set sf_xyz [measure center $sf]
set sf_xy  [lrange $sf_xyz 0 1]
set sf_z   [lindex $sf_xyz 2]
# now get the starting state which is identical to the old list
for {set i 0} {$i < $n} {incr i} {
  set ion [atomselect top "index [lindex $index $i]"]
  set xyz [measure center $ion]
  set z  [lindex $xyz 2]
  set xy [lrange $xyz 0 1] 
  if {$z > [expr $top - $dtop]} {
    lappend state 2; # its above
  } elseif {$z < [expr $bottom+$dbot]} {
    lappend state -2; # its below
  } else {
    if {[vecdist $sf_xy $xy] < $radius} {
      lappend state 0; # its in the channel
    } else {
      if {$z > $sf_z} {
        lappend state 2;  # its the xy plane but not in the channel. Its above
      } else {
        lappend state -2; # its the xy plane but not in the channel. Its below
      }
    }
  }
#puts $log "state $i:  [lindex $state 0]"
}


# the total event counter
set up 0
set down 0
# now go thru the frames
proc analysis {frame} {
  global out index resid n state k sf dt up down radius log dtop dbot
  # calculate the time
  set time [expr $frame * $dt]
  # the selectivity filter
  $sf frame $frame
  # its edges
  set top    [lindex [lindex [measure minmax $sf] 1] 2]
  set bottom [lindex [lindex [measure minmax $sf] 0] 2]
  set sf_xyz [measure center $sf]
  set sf_xy  [lrange $sf_xyz 0 1]
  set sf_z   [lindex $sf_xyz 2]
  # now get the positions of the ions
  $k frame $frame
  set i 0
  foreach xyz [$k get {x y z}] {
    # check to see if it near the selectivity filter
    set xy [lrange $xyz 0 1] 
    set z  [lindex $xyz 2]
#	puts $log "top: $top	bot: $bottom	XYDIST: [vecdist $sf_xyz $xyz]	sfz: $sf_z	z: $z"	
	
    if {$z > [expr $top-$dtop]} {
	   if {[lindex $state $i] == -1} {
          incr up
          puts $log "Event up Water index [lindex $index $i] resid [lindex $resid $i] at frame $frame"
          puts $out "$time [lindex $index $i] $up $down"
        }
        lset state $i 2
    } elseif {$z < [expr $bottom+$dbot]} {
	   if {[lindex $state $i] == 1} {
          incr down
          puts $log "Event down Water index [lindex $index $i] resid [lindex $resid $i] at frame $frame"
          puts $out "$time [lindex $index $i] $up $down"
        }
        lset state $i -2
      } else {
        if {[vecdist $sf_xy $xy] < $radius} {
	  		if	{[lindex $state $i] ==  2} {
            lset state $i 1
          } elseif {[lindex $state $i] ==  1} {
            lset state $i 1
          } elseif {[lindex $state $i] ==  0} {
            lset state $i 0
          } elseif {[lindex $state $i] == -1} {
            lset state $i -1
          } elseif {[lindex $state $i] == -2} {
            lset state $i -1
          }
        } else {
      	  if {$z > [expr $top-$dtop]} {	#$sf_z
		  	if	{[lindex $state $i] ==  -2} {
	      	    lset state $i 2;  # its the xy plane but not in the channel. Its above
			}
      	  } elseif {$z < [expr $bottom+$dbot]} {
		  	  if	{[lindex $state $i] ==  2} {
	      	    lset state $i -2; # its the xy plane but not in the channel. Its below
			}
      	  }
        }
      }
#	puts $log "frame  $frame - state $i = [lindex $state $i]"
    incr i
  }
  report $frame 100 "Down: $down Up: $up"
}
proc end {} {
  global out log
  close $out
  close $log
  exit
}
#waitfor [lindex $xtcs end] finish
#eval "bigdcd dcd analysis $xtcs"
waitfor $trj end
bigdcd $trj_type analysis $trj

