# this script is aimed at find water connectivity from a particualr posistion in the protein to 
# the bulk. The bulk water is defined by a cutoff in terms of z. If the protein is centered
# in 0 0 0 then the bluck water will be roughly at +/- 30A or which ever parameter you give.
# I ma only looking at water oxygens


# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/simple_water_file.tcl -args after_6000.pdb system.xtc result.txt 8 4 15 water_file.txt

source ~arkini/data/scripts/Tcl/bigdcd.tcl
source ~arkini/data/scripts/Tcl/waitfor.tcl
source ~arkini/data/vmd_scripts/shy_utilities.tcl

# the command line arguments parsing
set pdb    	[lindex $argv 0] ;#
set xtc    	[lindex $argv 1] ;#
set file   	[lindex $argv 2] ;#
set wsd         [lindex $argv 3]; # how close are water moleculaes to the site of interest
set wwd         [lindex $argv 4]; # how close are two water moleculaes to be defined as connected
set z_threshold [lindex $argv 5] ;# when the water file has reachedbulk water
set file        [lindex $argv 6] ;# file name

# load the protein  
mol load pdb $pdb  
# open the file handle
set fid [open $file w]

# run thru the trajectory
proc analyze {frame} {
  global z_threshold wsd wwd fid
  center_by_selection "protein"
  report $frame 50
  # first do 163 
  set sel "resid 155 and type OD1 OD2"
  # start of recursive selections
  set water_1  "type OWS and within $wsd of  ($sel)"
  set water_2  "type OWS and within $wwd of  ($water_1)"
  set water_3  "type OWS and within $wwd of  ($water_2)"
  set water_4  "type OWS and within $wwd of  ($water_3)"
  set water_5  "type OWS and within $wwd of  ($water_4)"
  set water_6  "type OWS and within $wwd of  ($water_5)"
  set water_7  "type OWS and within $wwd of  ($water_6)"
  set water_8  "type OWS and within $wwd of  ($water_7)"
  set water_9  "type OWS and within $wwd of  ($water_8)"
  set water_10 "type OWS and within $wwd of ($water_9)"
  set final_selection [atomselect top $water_9]
  # end of recursive selections
  set minmax [measure minmax $final_selection]
  set z_max [lindex [lindex $minmax 1] 2]
  set z_min [lindex [lindex $minmax 0] 2]
  if {$z_max  >= $z_threshold} {
    puts $fid "Asp 163 in frame $frame is open to cytoplasm"
    puts "Asp 163 in frame $frame is open to cytoplasm"
  }
  if {[expr 0 - $z_threshold] >= $z_min} {
    puts $fid "Asp 163 in frame $frame is open to periplasm"
    puts "Asp 163 in frame $frame is open to cytoplasm"
  }
  # now do 164
  set sel "resid 156 and type OD1 OD2"
  # start of recursive selections
  set water_1  "type OWS and within $wsd of  ($sel)"
  set water_2  "type OWS and within $wwd of  ($water_1)"
  set water_3  "type OWS and within $wwd of  ($water_2)"
  set water_4  "type OWS and within $wwd of  ($water_3)"
  set water_5  "type OWS and within $wwd of  ($water_4)"
  set water_6  "type OWS and within $wwd of  ($water_5)"
  set water_7  "type OWS and within $wwd of  ($water_6)"
  set water_8  "type OWS and within $wwd of  ($water_7)"
  set water_9  "type OWS and within $wwd of  ($water_8)"
  set water_10 "type OWS and within $wwd of ($water_9)"
  set final_selection [atomselect top $water_9]
  # end of recursive selections
  set minmax [measure minmax $final_selection]
  set z_max [lindex [lindex $minmax 1] 2]
  set z_min [lindex [lindex $minmax 0] 2]
  if {$z_max  >= $z_threshold} {
    puts $fid "Asp 164 in frame $frame is open to cytoplasm"
    puts "Asp 164 in frame $frame is open to cytoplasm"
  }
  if {[expr 0 - $z_threshold] >= $z_min} {
    puts $fid "Asp 164 in frame $frame is open to periplasm"
    puts "Asp 164 in frame $frame is open to cytoplasm"
  }
}

proc end {} {
  exit
}
waitfor $xtc end
bigdcd xtc analyze $xtc















