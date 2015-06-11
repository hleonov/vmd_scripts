# this script should run through a trajectory and calculate the z posion of the ion relative to its initialposition (taken as 0

# it runs on a single xtc and takes a while since it uses bigdcd to read one frame at a time

# run this script as follows:
# vmd -dispdev text -e ~/data/vmd_scripts/ion_z_relative_to_strat.tcl -args system_from_em.pdb system.xtc results.xls Na

if {[llength $argv] < 3} {
  puts "Usage example:"
  puts "vmd -dispdev text -e ~/data/vmd_scripts/ion_z_relative_to_strat.tcl -args system_from_em.pdb system.xtc results.xls Na"
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
set ion [atomselect top $target frame 0]
# computations at frame 0
set ion_at_0 [lindex [measure center $ion] 2]
set ref [atomselect top "protein and name CA" frame 0]

proc z_diff {frame} {
  global ion ref ion_at_0 fid
  # update frame
  $ion frame $frame
  set all [atomselect top "all"]
  set sel [atomselect top "protein and name CA"]
  #$all move [measure fit $sel $ref]
  
  set ion_now [lindex [measure center $ion] 2]
  set z_diff [expr $ion_now - $ion_at_0]  
  # output
  puts $fid "$frame $z_diff"
  # this makes a nice printout every 25 frames
  if {[expr $frame / 25.0] == int([expr $frame / 25.0])} {
   puts "Analyzing frame: $frame => $z_diff"
  }
}
proc end {} {
  global fid
  close $fid
  exit
}
waitfor $xtc end
bigdcd xtc z_diff $xtc



















