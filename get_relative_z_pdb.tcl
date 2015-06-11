# this script should open a pdb and write out the relative
# position of the ion (input) and the lipid bilayer

# vmd -dispdev text -e ~/data/vmd_scripts/get_relative_z_pdb.tcl -args after_24000.pdb results.dat Na

if {[llength $argv] < 3} {
  puts "Usage example:"
  puts "vmd -dispdev text -e ~/data/vmd_scripts/get_relative_z_pdb.tcl -args after_24000.pdb results.dat Na"
  exit
}
# the command line arguments parsing
set pdb      [lindex $argv 0] ;
set file     [lindex $argv 1] ;
set target   [lindex $argv 2] ;
# set the target
set target "resname $target"
# here I load the protein  
mol load pdb $pdb  
# open the file handle
set fid [open $file w]
# selections
set ion [atomselect top $target]
set lipid [atomselect top "resname POP"]
# computation
set d [expr [lindex [measure center $ion weight mass] 2] - [lindex [measure center $lipid weight mass] 2] ]

puts $fid $d
puts  $d

exit















