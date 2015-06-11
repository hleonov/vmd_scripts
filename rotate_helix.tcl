source ~stan/Bin/utilities.tcl
set number_of_amino_acids_to_average 4
set buffer_from_end 1


proc rotate_helix {angle} {
  global number_of_amino_acids_to_average buffer_from_end
  set all [atomselect top all];
  set CAs [atomselect top "name CA"];
  # find the range of residues 
  set first_resid [lindex [$CAs get {resid}] 0];
  set last_resid  [lindex [$CAs get {resid}] end];
  # calculate the helix axis
  set a [expr $first_resid + $buffer_from_end]
  set b [expr $first_resid + $buffer_from_end + $number_of_amino_acids_to_average]
  set start [atomselect top "name CA and resid $a to $b"];
  set a [expr $last_resid - $buffer_from_end - $number_of_amino_acids_to_average]
  set b [expr $last_resid - $buffer_from_end]
  set end  [atomselect top "name CA and resid $a to $b"];
  set start_coords [measure center $start];
  set end_coords   [measure center $end  ];
  set matrix [trans bond $start_coords $end_coords $angle deg]
  $all move $matrix
}










