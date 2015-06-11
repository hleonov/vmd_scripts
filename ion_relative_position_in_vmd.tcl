proc z_diff { ion } {
  set sel [atomselect top "resname $ion" ]
  set lipid [atomselect top "resname POP"]
  set ion_lipid [expr [lindex [measure center $sel] 2] - [lindex [measure center $lipid] 2] ]
  set sel [atomselect top "resname $ion" frame 0]
  set lipid [atomselect top "resname POP" frame 0 ]
  set ion_lipid_0 [expr [lindex [measure center $sel] 2] - [lindex [measure center $lipid] 2] ]
  set diff [expr $ion_lipid - $ion_lipid_0]
  puts "At 0: $ion_lipid_0; Now: $ion_lipid; Difference: $diff"
}
