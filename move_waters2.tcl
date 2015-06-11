# this script will take 3 names~:
# 1 the file with the HOH whcih has only oxygens
# 2 a file with a single water now with hydrogens
# 3 an output file name
# it will move the waters to the location of the oxygens and write it out.

source /Users/hleonov/vmd_scripts/utilities.tcl
check_for_help "Usage: vmdt -e move_waters2.tcl -args -rf system_hole.pdb sol.pdb hoh_moved.pdb"
if {[lindex $argv 0] == "-rf"} {
  set pdb  [lindex $argv 1];
  set sol  [lindex $argv 2];
  set out  [lindex $argv 3];

} else {
  set pdb  [lindex $argv 0];
  set sol  [lindex $argv 1];
  set out  [lindex $argv 2];
}
# open the file handle
#set out [open [out_file $out] w];# check if the output file exists
# here I load the hoh water  
mol load pdb $pdb  
# now load the single water
mol load pdb $sol


set hoh_sel [atomselect 0 "resname HOH"]
set sol_sel [atomselect 1 "resname SOL"]
puts [$sol_sel num]

set i 0
set r 3784
foreach xyz [$hoh_sel get {x y z}] {
  $sol_sel moveby [vecinvert [vecsub [measure center $sol_sel] $xyz]]
  $sol_sel set resid $r
  $sol_sel set resname HOH
  $sol_sel writepdb hoh_moved_${i}.pdb
  incr i 
  incr r
}

set cmd "cat hoh_moved_*.pdb | grep ATOM > $out"
exec sh -c $cmd 

exit

