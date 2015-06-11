# this script should run through a trajectory and generate a histogram with waters
# the outcome is waters per A cubed. Bulck water density is 0.033.
set pdb_file_name closewaters.pdb
set residue 156
set file_name results_164.xls
set number_of_dcds 17
set dcd_name "dcds/close"
set range 35.0
# here I load the protein  
mol load pdb $pdb_file_name  
# define some basic parameters
set framestotal 0
set number_of_slices 50.0
set radius 3
set total_number_frames 0
# some simple calculations
set slab_width [expr 2 * $range * 1.0 / $number_of_slices ]
set r2 [expr $radius*$radius]
# open the file handle
set fid [open $file_name w]
# initialze the counting array
for {set i 0} {$i <= $number_of_slices} {incr i} {
  set counting_array($i) 0
}
# read each dcd on its own.
for {set j 0} {$j < $number_of_dcds} {incr j} {
  set dcd_file "$dcd_name$j.dcd"
  animate read dcd $dcd_file waitfor all
  # the analysis
  # run thru the frames
  set number_of_frames [molinfo top get numframes]
  # Center around the selection
  set start [atomselect top "all" frame 0]
  set sel [atomselect top "resid $residue and protein and type OD1 OD2" frame 0]
  $start moveby [vecinv [measure center $sel]]
  # Move it such that the minimal z coordiante (of future waters) is at 0
  $start moveby [list 0 0 $range]
  # ste the reference for further alignment
  set ref [atomselect top "protein and name CA" frame 0]
  for {set k 0} {$k < $number_of_frames} {incr k} { 
    # Align and move
    set all [atomselect top "all" frame $k]
    set ca [atomselect top "protein and name CA" frame $k]
    $all move [measure fit $ca $ref]
    # Now define the selections for the water location
    set water_sel [atomselect top "name OWS and x^2 + y^2 <= $r2 and (z >= 0 ) and (z <= 2* $range )" frame $k]
    # now count the number of atoms in each histogram bin
    foreach z [$water_sel get {z}] {     
      set slice [ expr floor( [expr $z / $slab_width ] ) + 1 ] 
      set slice [ expr int($slice)] 
      incr counting_array($slice)
    }
    incr total_number_frames
  } 
  # delete the trajectory file before loading the next one (but not the initial frame which is the PDB file)
  animate delete beg 1
}
# calculate the slab volume
set slab_volume [expr $r2 * 3.14159265358979 * $slab_width]
puts $total_number_frames
for {set l 1} {$l < $number_of_slices} {incr l} { 
  set counting_array($l) [expr $counting_array($l) * 1.0 / ( $total_number_frames * $slab_volume ) ]
  puts $fid "$l $counting_array($l)"
}
close $fid
exit

