# this doesn;t work yet!!!!!


# this script is aimed at find water connectivity from a particualr posistion in the protein to 
# the bulk. The bulk water is defined by a cutoff in terms of z. If the protein is centered
# in 0 0 0 then the bluck water will be roughly at +/- 30A or which ever parameter you give.
# I ma only looking at water oxygens


# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/water_file.tcl -args after_6000.pdb system.xtc result.xls 156

source ~arkini/data/scripts/Tcl/bigdcd.tcl
source ~arkini/data/scripts/Tcl/waitfor.tcl
source ~arkini/data/vmd_scripts/shy_utilities.tcl

# the command line arguments parsing
set pdb   [lindex $argv 0] ;# 
set xtc   [lindex $argv 1] ;# 
set file  [lindex $argv 2] ;# 
set resid [lindex $argv 3]; # the resid to start searchign from


# some basic parameters 
set closeness   12;  # how close are two water moleculaes to be defined as connected
set z_threshold 15; # when the water file has reachedbulk water
# load the protein  
mol load pdb $pdb  

proc find_close_waters {sel frame} {
  global closeness
  set near [atomselect top "type OWS and within $closeness of [$sel text]"]
  return [$near list]
}

# run thru the trajectory
proc analyze {frame} {
  global resid z_threshold
  # center the protein at 0 0 0
  center_by_selection "protein"
  # puts "====================================================="
  # puts "Frame => $frame"
  # set the selection to what I want
  set sel [atomselect top "resid $resid and type OD1 OD2"]
  # define the global water list (this isn't a selection it is a list)
  # and in it put the selection atoms although they aren't water
  set global_water_list [$sel list]
  # report every nth frames
  report $frame 50
  # generate the initial list fo close waters
  # set the selection to what I want
  set sel [atomselect top "index [$sel list]"]
  # find the close water and put them into a list
  set close_waters [find_close_waters $sel $frame]
  set counter 0
  # now searching a continuous water file
  while {[llength $close_waters] > 0} {
    incr counter
    #puts "$counter global_water_list: $global_water_list"
    #puts "$counter close_waters: $close_waters"
    set sel [atomselect top "index [$sel list]"]
    # find the close waters and put them into a list
    set close_waters [find_close_waters $sel $frame]
    # now check to see if any of these new waters are new (i.e. not present in $global_water_list)
    foreach water $close_waters {
      # this assumes that the list has no duplications, which it shouldn't
      set i [lsearch $close_waters $water]
      foreach global_water $global_water_list {
  	if {$water == $global_water} {
  	  set close_waters [lreplace $close_waters $i $i]
  	}
      }
    }
    # check to see if the list has any element thta are more than +/- z form the center of the protein
    set breakout 0
    foreach j $close_waters {
      set tmp_sel [atomselect top "index $j"]
      set z [$tmp_sel get {z}]
      if {$z >= $z_threshold || [expr 0 - $z_threshold] >= $z } {
        set breakout 1
      }
    } 
    if {$breakout > 0} {
      puts "Frame $frame has a contiuous file"
      break
    }   
    # check to see if the list is empty
    # if it is then we have failed to fine a water file
    if {[llength $close_waters] == 0} {
      # puts "No water file at frame $frame"
    }
    # add the none zero list of close waters to the global list of waters
    foreach i $close_waters {
      lappend global_water_list $i
    }
  }
  
}

proc end {} {
  exit
}
waitfor $xtc end
bigdcd xtc analyze $xtc















