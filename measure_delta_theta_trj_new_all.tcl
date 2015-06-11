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


proc tdm_atoms {resid frame} {
  # assuem that the tdm attom are the clostst to the C
  set tdm1 [atomselect top "name TDM1 and within 1.5 of (resid $resid and type C)" frame $frame]
  set tdm2 [atomselect top "name TDM2 and within 1.5 of (resid $resid and type C)" frame $frame]
  if {[$tdm1 num] < 1} {
    return "-"
  } else {
    set tdm1_coords [measure center $tdm1]
    set tdm2_coords [measure center $tdm2]
    set vector [vecsub $tdm1_coords $tdm2_coords]
    set z [list 0 0 1]
    set angle [angle_between_two_vectors $vector {0 0 1}]
    if {$angle > 90} {
      return [expr 180 - $angle]
    } else {
      return $angle
    }
  }
}
proc tdm_real {resid frame} {
  set c [atomselect top "resid $resid and type C" frame $frame]
  set o [atomselect top "resid $resid and type O" frame $frame]
  set n [atomselect top "resid [expr $resid + 1] and type N" frame $frame]
  set c_coords [measure center $c]
  set o_coords [measure center $o]
  set n_coords [measure center $n]
  set c_o [vecnorm [vecsub $o_coords $c_coords]];
  set tdm1 [vecadd $c_coords [vecscale $c_o 0.868]]
  set c_n [vecnorm [vecsub $n_coords $c_coords]];
  set tdm2 [vecadd $c_coords [vecscale $c_n 0.505]]
  set vector [vecsub $tdm1 $tdm2]
  set z [list 0 0 1]
  set angle [angle_between_two_vectors $vector {0 0 1}]
  if {$angle > 90} {
    return [expr 180 - $angle]
  } else {
    return $angle
  }
}
# Now run it.

set residues {5 6 7 8 9 12 15 17 19 22}
set true {16.5 30.985 52.52 49.35 20.95 33.095 54.875 47.685 47.35 63.865}; # low pH
# set true {51.18 49.97 53.83 50.95 51.88 57.14 55.97 44.65 52.42 56.82}; # high pH
# set true {20.689 39.895 69.357 53.125 4.1486 36.573 62.77 41.461 49.28 68.43}; # NMR

foreach res_num $residues {
  # puts -nonewline "\t $res_num \t $res_num"
#  puts -nonewline $out "\t $res_num \t $res_num"
  puts -nonewline $out "\t $res_num"
}
puts $out " "
for {set i 0} {$i <= [molinfo top get numframes]} {incr i} {
  puts -nonewline $out "$i"
  # puts -nonewline  "$i"
  set j 0;
  foreach res_num $residues {
     set current [lindex $true $j]
#     set a [expr abs([expr [tdm_atoms $res_num $i] - $current])]
     set b [expr abs([expr [tdm_real $res_num $i] - $current])]
#     puts -nonewline $out "\t $a\t $b" 
     puts -nonewline $out "\t $b"   
     set j [expr $j + 1]
  #  puts -nonewline   "\t [tdm_atoms $res_num $i]\t [tdm_real $res_num $i]"
  }
  # puts      " "
  puts $out " "
}
puts " "
puts "Finished"
close $out
exit









