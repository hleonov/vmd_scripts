# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/water_shell.tcl -args after_24000.pdb system.xtc water_shell.xls K 4
source ~arkini/data/scripts/Tcl/bigdcd.tcl
source ~arkini/data/scripts/Tcl/waitfor.tcl
source ~arkini/data/vmd_scripts/shy_utilities.tcl

# the command line arguments parsing
set pdb      [lindex $argv 0] ;
set xtc      [lindex $argv 1] ;
set file     [lindex $argv 2] ;
set ion      [lindex $argv 3] ;
set distance [lindex $argv 4] ;
# here I load the protein  
mol load pdb $pdb
# open the file handle
set fid [open $file w]
# run thru the trajectory
proc analyze {frame} {
  global sel file ion distance fid
  report $frame 100
  set sel [atomselect top "resname $ion"]
  set water [atomselect top "same residue as type OWS and within $distance of resname $ion"]
  puts $fid [expr [$water num] / 3]
}
proc end {} {
  exit
}
waitfor $xtc end
bigdcd xtc analyze $xtc















