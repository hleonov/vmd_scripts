proc write_selective_trj {mol} {
  set number_of_frames [expr [molinfo $mol get numframes] - 1]
  set long_list {}
  set sel [atomselect $mol "(same residue as water and within 4 of protein) or (same residue as resname POP and within 4 of protein)"]
  for {set i 0} {$i <= $number_of_frames} {incr i} {
    $sel update
    $sel frame $i
    foreach elem [$sel list] {lappend long_list $elem}
    set long_list [lsort -uniq $long_list]
  }
  set sel_text "index $long_list"
  set sel [atomselect $mol "protein or $sel_text"]
  puts $number_of_frames
  animate write trr selective_system.trr beg 0 end $number_of_frames sel $sel $mol
  $sel writepdb selective_system.pdb
}
