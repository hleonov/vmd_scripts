source shy_utilities.tcl
package require autoionize


# the -rf needs to be the first argument due to the ability to grab many xtcs
check_for_help "Usage: vmdt -e ~/TclScripts/surface_area_per_lipid.tcl -args -rf system_pr_500.pdb results.xls trjformat MCnumber system_pr_0.trr
trjformat is the type of the trajectory file format (e.g. xtc).
MCnumber is the number of trials to ue in the Monte Carlo
Example:
vmdt -e ~/TclScripts/surface_area_per_lipid.tcl -args -rf system_pr_500.pdb results.xls trr 10 system_pr_0.trr"
if {[lindex $argv 0] == "-rf"} {
  set pdb  [lindex $argv 1];
  set out  [lindex $argv 2];
  set TrjName  [lindex $argv 3];
  set NumberOfTrials  [lindex $argv 4];
  set xtcs [lrange $argv 5 end]

} else {
  set pdb  [lindex $argv 0];
  set out  [lindex $argv 1];
  set TrjName  [lindex $argv 2];
  set NumberOfTrials  [lindex $argv 3];
  set xtcs [lrange $argv 4 end]
}
# here I load the protein  
mol load gro $pdb  
# open the file handle
set out [open [out_file $out] w]; # see -rf flag note above

# Defenitions
set dt 10; # the time of each frame in ps
set AtomName "P"; # the name of the lipid atom used to get the xy plane
set d 2; # distance from protein to check

# pick the last atom in the system. It is normamly an ion or part of water
set sel [atomselect top "all"]
set last [lindex [$sel list] end]
set sel [atomselect top "index $last"]
# more selections
set water [atomselect top "water"]


proc surface_area {} {
  global AtomName water NumberOfTrials sel d last out
  # move the system such that the system center is at 0 0 0
  center_by_selection "all"
  # measure the average z of all the phosphates of the upper leafleat
  set top    [atomselect top "type $AtomName and z > 0"]
  set bottom [atomselect top "type $AtomName and z <= 0"]
  # the number of lipids at tope and at bottom
  set lipids_at_top    [$top num]
  set lipids_at_bottom [$bottom num]
  # figure out the average z of each leaflet
  set z_top    [lindex [measure center $top]    2]
  set z_bottom [lindex [measure center $bottom] 2]
  # puts "[$top num] Lipids in the upper leaflet (z = $z_top)"
  # puts "[$bottom num] Lipids in the upper leaflet (z = $z_bottom)"
  # now figure out the range that x and Y can take
  # this is the size of the waters
  set xmin [lindex [lindex [measure minmax $water] 0] 0]
  set ymin [lindex [lindex [measure minmax $water] 0] 1]
  set xmax [lindex [lindex [measure minmax $water] 1] 0]
  set ymax [lindex [lindex [measure minmax $water] 1] 1]
  # puts "Xmin: $xmin; Xmax: $xmax; Ymin: $ymin; Ymax:$ymax"
  set xrange [expr $xmax - $xmin]
  set yrange [expr $ymax - $ymin]
  set area [expr $xrange * $yrange]
  # puts "$xrange $yrange"
  set top 0
  set bottom 0
  for {set i 1} {$i <= $NumberOfTrials} {incr i} {
    set xrnd [expr rand() * $xrange];
    set yrnd [expr rand() * $xrange];
    set x [expr $xmax - $xrnd];
    set y [expr $ymax - $yrnd];
    # the top layer
    $sel moveto [list $x $y $z_top]
    set close_to_protein [atomselect top "protein and within $d of index $last"]
    if {[$close_to_protein num] > 0} {
      incr top
    }    
    # the bottom layer
    $sel moveto [list $x $y $z_bottom]
    set close_to_protein [atomselect top "protein and within $d of index $last"]
    if {[$close_to_protein num] > 0} {
      incr bottom
    }    
  }
  # the following are the ratio fo the lipids within the membrane
  set top    [expr 1 - ($top * 1.0    / $NumberOfTrials)]
  set bottom [expr 1 - ($bottom * 1.0 / $NumberOfTrials)]
  set area_per_lipid_top    [expr $top    * $area / $lipids_at_top]
  set area_per_lipid_bottom [expr $bottom * $area / $lipids_at_bottom]
  #puts $out "$area_per_lipid_top\t $area_per_lipid_bottom"
  #puts "Top: $area_per_lipid_top; Bottom: $area_per_lipid_bottom"
  
  return [list $area_per_lipid_top $area_per_lipid_bottom]

}

# now go thru the frames
proc analysis {frame} {
  global AtomName water NumberOfTrials sel d last out dt
  # calculate the time
  set time [expr $frame * $dt]
  set top_bottom [surface_area]
  set top [lindex $top_bottom 0]
  set bottom [lindex $top_bottom 1]
  puts $out "$frame $top $bottom"
  report $frame 10 "$time $top $bottom"
}
proc finish {} {
  global out
  close $out
  exit
}
waitfor [lindex $xtcs end] finish
eval "bigdcd $TrjName analysis $xtcs"
