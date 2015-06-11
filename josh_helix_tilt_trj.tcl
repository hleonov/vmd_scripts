# this will read in a a pdb and will output
source ~stan/Bin/utilities.tcl

if {[lindex $argv 0] == "-rf"} {
  set out  [lindex $argv 1];
  set pdb  [lindex $argv 2];
  set xtc  [lindex $argv 3];
} else {
  set out  [lindex $argv 0];
  set pdb  [lindex $argv 1];
  set xtc  [lindex $argv 2];
}
mol load pdb $pdb;
animate read xtc $xtc waitfor all;
set out  [open $out w];



proc tilt {frame} {
  set first_resid [find_first_residue];
  set last_resid  [find_end_residue];
  set start [atomselect top "resid $first_resid to [expr $first_resid + 3] and name CA" frame $frame]
  set end   [atomselect top "resid $last_resid to[expr $last_resid  - 3] and name CA" frame $frame]
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
proc find_first_residue { } {
  set CAs [atomselect top "name CA"];
  # find the range of residues 
  set first_resid [lindex [$CAs get {resid}] 0];
  return $first_resid
}
proc find_end_residue { } {
  set CAs [atomselect top "name CA"];
  set last_resid  [lindex [$CAs get {resid}] end];
  return $last_resid
}

for {set i 0} {$i <= [molinfo top get numframes]} {incr i} {
  puts -nonewline $out "$i\t"
  puts $out [tilt $i]

  #puts -nonewline "$i\t"
  #puts  [tilt $i]
}
close $out
exit





