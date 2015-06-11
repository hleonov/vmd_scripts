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
set closeness   10;  # how close are two water moleculaes to be defined as connected
set z_threshold 15; # when the water file has reachedbulk water





# load the protein  
mol load pdb $pdb  

# define the global water list (this isn't a selection it is a list)
global global_water_array
foreach k {39253 41269} {
  set global_water_array($k) 1
}
  

proc find_close_waters {sel frame} {
  global closeness
  set near [atomselect top "type OWS and within $closeness of [$sel text]"]
  return [$near list]
}


# run thru the trajectory
proc analyze {frame} {
  global resid global_water_array
  #report $frame 50
  # set the selection to what I want
  set sel [atomselect top "resid $resid and protein"]
  # find the close water and put them into a new list
  set close_waters [find_close_waters $sel $frame]
  
  foreach close_water $close_waters {
  	set close_water_array($close_water) 1
  }
  
  #puts [llength $close_waters]
  # now check to see if any of these new waters are new (i.e. not present in $water_list)
  # if they are new add it to the list

  foreach close_water $close_waters {
    if {[info exists global_water_array($close_water)]} {
	# using a list
	set i [lsearch $close_waters $close_water]
	set close_waters [lreplace $close_waters $i $i]
	#  using an array
	#unset close_water_array($close_water)
    }
  }
  puts $close_waters
  #set close_waters [array names close_water_array]
  
}

proc end {} {
  exit
}
waitfor $xtc end
bigdcd xtc analyze $xtc















