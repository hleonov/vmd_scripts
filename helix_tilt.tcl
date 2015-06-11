# this will read in a a pdb and will output
# the tilt of the helix

source ~/vmd_scripts/utilities.tcl
#check_for_help "Usage: vmdt -e ~/Bin/helix_tilt.tcl -args my_proein.pdb"

set pdb    [lindex $argv 0];
set resb   [lindex $argv 1];
set rese   [lindex $argv 2];

mol load pdb $pdb  


proc tilt {b e} {

#  set first_resid [lindex [find_end_residues] 0];
#  set last_resid  [lindex [find_end_residues] end];
	set first_resid $b
	set last_resid $e
  set first_resid [expr $first_resid + 3]
  set last_resid  [expr $last_resid - 3]
  # set start [atomselect top "resid $first_resid [expr $first_resid + 4] and name CA"]
  set start [atomselect top "resid $b to $first_resid and name CA"]
  set end   [atomselect top "resid $last_resid  to $e and name CA"]
  set start_coords [measure center $start]
  set end_coords   [measure center $end]
  set helix [vecsub $start_coords $end_coords]
  set angle [angle_between_two_vectors $helix {0 0 1}]
  if { $angle > 90} {
    return [expr 180 - $angle]
  } else {
    return $angle
  }
}
proc find_end_residues { } {
  set CAs [atomselect top "type CA"];
  # find the range of residues 
  set first_resid [lindex [$CAs get {resid}] 0];
  set last_resid  [lindex [$CAs get {resid}] end];
  return [list $first_resid $last_resid]
}

set a2  [angle_between_two_vectors {1 0 0} {0 0 1}]
puts "90 degrees: $a2" 
puts " ----------------- Tilt angle ------------------ "
puts [tilt $resb $rese]
puts " ----------------------------------------------- "
exit


