source ~stan/Bin/utilities.tcl
check_for_help "Usage: vmdt -e ~/Bin/place_princ.tcl -args -rf my_proein.pdb output_file.txt"
if {[lindex $argv 0] == "-rf"} {
  set pdb    [lindex $argv 1];
  set out    [lindex $argv 2];
} else {
  set pdb    [lindex $argv 0];
  set out    [lindex $argv 1];
}
mol load pdb $pdb
set out [open [out $out] w];


proc distance {resid_aa resid_prnc chain} {
  if {$chain > 0} {
    set ca1  [atomselect top "resid $resid_aa and chain $chain and type CA"];
    set ca2  [atomselect top "resid [expr $resid_aa + 1] and chain $chain and type CA"];
    set ca3  [atomselect top "resid [expr $resid_aa + 2] and chain $chain and type CA"];
    set ca4  [atomselect top "resid [expr $resid_aa + 3] and chain $chain and type CA"];
    set prnc [atomselect top "resid $resid_prnc and chain $chain and type PRNC"];
    
  } else {
    set ca1  [atomselect top "resid $resid_aa and type CA"];
    set ca2  [atomselect top "resid [expr $resid_aa + 1] and type CA"];
    set ca3  [atomselect top "resid [expr $resid_aa + 2] and type CA"];
    set ca4  [atomselect top "resid [expr $resid_aa + 3] and type CA"];
    set prnc [atomselect top "resid $resid_prnc and type PRNC"];
  }
  set dis1 [vecdist [measure center $prnc] [measure center $ca1]]
  set dis2 [vecdist [measure center $prnc] [measure center $ca2]]
  set dis3 [vecdist [measure center $prnc] [measure center $ca3]]
  set dis4 [vecdist [measure center $prnc] [measure center $ca4]]
  
  return [list $dis1 $dis2 $dis3 $dis4]
}

# Aux procedures
proc find_end_residues { } {
  set CAs [atomselect top "all and type CA"];
  # find the range of residues 
  set first_resid [lindex [$CAs get {resid}] 0];
  set last_resid  [lindex [$CAs get {resid}] end];
  return [list $first_resid $last_resid];
}

proc find_very_end_residues { } {
  set all [atomselect top "type TDM2"];
  # find the range of residues 
  set first_resid [lindex [$all get {resid}] 0];
  set last_resid  [lindex [$all get {resid}] end];
  return [list $first_resid $last_resid];
}

proc find_end_id { } {
  # remeber the index is zero base
  set all [atomselect top "all"];
  set number_of_atoms [$all num];
  return [expr $number_of_atoms - 1]
}

proc sum {args} {
  set res 0
  foreach n $args {
    set res [expr $res + $n];
  }
  return $res 
}

proc average {args} {
 puts $args
 return [expr [eval sum $args] / [llength $args]]
}

#running it

set first_res [lindex [find_end_residues] 0];
set end_res   [lindex [find_end_residues] 1];
puts [find_very_end_residues]
set end_all   [lindex [find_very_end_residues] 1];
# Taking all the way until the last 4, and -1 to exclude the last one without a partner
set end [expr $end_res - 4];


for {set i $first_res} {$i < $end} {incr i} {
   set exp [distance $i [expr $i + $end_all] 0];
   set dis1 [lindex $exp 0]
   set dis2 [lindex $exp 1]
   set dis3 [lindex $exp 2]
   set dis4 [lindex $exp 3]
   # puts "diss: $dis1 $dis2 $dis3 $dis4"
   lappend d1 $dis1
   lappend d2 $dis2
   lappend d3 $dis3
   lappend d4 $dis4
   # puts $out "$exp"
}

set ave1 [expr [eval sum $d1] / [llength $d1]]
set ave2 [expr [eval sum $d2] / [llength $d2]]
set ave3 [expr [eval sum $d3] / [llength $d3]]
set ave4 [expr [eval sum $d4] / [llength $d4]]

# set ave1 [average $d1]
# puts "ave1 is: $ave1"
# set ave2 [average $d2]
# set ave3 [average $d3]
# aet ave4 [average $d4]


set exp2 [list $ave1 $ave2 $ave3 $ave4]
puts $out "$exp2"

puts "completed"


exit



  

