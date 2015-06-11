# run as: vmd -dispdev text -e fold.tcl
source ~/shy_utilities.tcl
set in  [open "data/cath_35_list_swiss_prot_location.txt_old" r];

# define the list for later on making a histogram
set data_cytoplasmic {}
set data_secreted    {}
while {[gets $in line] > -1} {
  set entry "data/CathDomainPdb.S35.v3.2.0/${line}"
  set pdb [lindex $entry 0]
  set location [lindex $entry 2]
  puts stderr "$pdb $location"
  mol delete all;
  if {$location == "Cytoplasm"} {
    mol load pdb $pdb;
    # now run thru every residue that is both not polar
    set sel [atomselect top "protein and (not polar) and name CA"]
    set protein [atomselect top "protein"]
    set buried {}
    foreach i [$sel list] resid [$sel get {resid}] {
      set residue [atomselect top "(same residue as index $i) and sidechain"]
      set area [area_per_resid $residue $protein 1.4]
      # calculate the accesibility ratio
      set resname  [[atomselect top "resid $resid and name CA"] get {resname}]
      set max [get_maximal_surface_area_of_residue $resname]
      if {$max < 1000} {
        set ratio [expr $area / $max]
        if {$ratio <= 0.1} {
          #puts "$i $resid $resname $area $max $ratio"
          lappend buried $resid
        }
      }
    }
    # now run thru the residues and print out the resid differences
    for {set i 1} {$i < [llength $buried]} {incr i} {
      set j [expr $i - 1]
      set d [expr [lindex $buried $i] - [lindex $buried $j]]
      lappend data_cytoplasmic $d
    }
  } elseif {$location == "Secreted"} {
    mol load pdb $pdb;
    # now run thru every residue that is both not polar
    set sel [atomselect top "protein and (not polar) and name CA"]
    set protein [atomselect top "protein"]
    set buried {}
    foreach i [$sel list] resid [$sel get {resid}] {
      set residue [atomselect top "(same residue as index $i) and sidechain"]
      set area [area_per_resid $residue $protein 1.4]
      # calculate the accesibility ratio
      set resname  [[atomselect top "resid $resid and name CA"] get {resname}]
      set max [get_maximal_surface_area_of_residue $resname]
      if {$max < 1000} {
        set ratio [expr $area / $max]
        if {$ratio <= 0.1} {
          #puts "$i $resid $resname $area $max $ratio"
          lappend buried $resid
        }
      }
    }
    # now run thru the residues and print out the resid differences
    for {set i 1} {$i < [llength $buried]} {incr i} {
      set j [expr $i - 1]
      set d [expr [lindex $buried $i] - [lindex $buried $j]]
      lappend data_secreted $d
    }
  }
}
print_x_y_lists [histogram $data_secreted 5] "result_secreted.dat"
print_x_y_lists [histogram $data_cytoplasmic 5] "result_cytoplasmic.dat"
print_x_y_lists [histogram $data_secreted 5] 0
print_x_y_lists [histogram $data_cytoplasmic 5] 0
exit;










