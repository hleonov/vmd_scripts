# this script should run through a trajectory and calculate the z posion of the Na

# it runs on a single xtc and takes a while since it uses bigdcd to read one frame at a time

# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/ion_relative_position.tcl -args after_24000.pdb system.xtc results.dat Na

if {[llength $argv] < 3} {
  puts "Usage example:"
  puts "vmd -dispdev text -e ~/data/vmd_scripts/ion_relative_position.tcl -args after_24000.pdb system.xtc results.dat Na"
  exit
}


source ~arkini/data/scripts/Tcl/bigdcd.tcl
source ~arkini/data/scripts/Tcl/waitfor.tcl

# the command line arguments parsing
set pdb      [lindex $argv 0] ;# start_na.pdb
set xtc      [lindex $argv 1] ;# all.xtc
set file     [lindex $argv 2] ;# results.dat
set target   [lindex $argv 3] ;# Na
# set the target
set target "resname $target"





# here I load the protein  
mol load pdb $pdb  
# open the file handle
set fid [open $file w]
# selections
set na [atomselect top $target frame 0]
set lipid [atomselect top "resname POP" frame 0]
set protein [atomselect top "protein" frame 0]
# computations at frame 0
set na_lipid_0 [expr [lindex [measure center $na] 2] - [lindex [measure center $lipid] 2] ]
set protein_lipid_0 [expr [lindex [measure center $protein] 2] - [lindex [measure center $lipid] 2] ]
set na_protein_0 [expr [lindex [measure center $na] 2] - [lindex [measure center $protein] 2] ]
proc z_diff {frame} {
  global target fid na lipid protein na_lipid_0 protein_lipid_0 na_protein_0
  # update frame
  $na frame $frame
  $lipid frame $frame
  $protein frame $frame
  # calculation Na - lipid
  set na_lipid [expr [lindex [measure center $na] 2] - [lindex [measure center $lipid] 2] ]
  set na_lipid [expr $na_lipid - $na_lipid_0]
  # calculation Na - protein
  set na_protein [expr [lindex [measure center $na] 2] - [lindex [measure center $protein] 2] ]
  #set na_protein [expr $na_protein - $na_protein_0]
  # calculation protein - lipid
  set protein_lipid [expr [lindex [measure center $protein] 2] - [lindex [measure center $lipid] 2] ]
  set protein_lipid [expr $protein_lipid - $protein_lipid_0]
  # output
  puts $fid "$frame $na_protein"
  # this makes a nice printout every 25 frames
  if {[expr $frame / 25.0] == int([expr $frame / 25.0])} {
   puts "Analyzing frame: $frame => $na_protein"
  }
}
proc end {} {
  global fid
  close $fid
  exit
}
waitfor $xtc end
bigdcd xtc z_diff $xtc



















