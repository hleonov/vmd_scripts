# run this from within vmd tkconsoule
# we use this to read the results of the h-bond analysis
# make sure that the skip number was the same

#cd /bigdisk2/Raphael/Li/all_charged/
#mol delete all
# mol load pdb md_0_8.pdb
# animate read xtc md_0_8.xtc skip 20 waitfor all
# source /Users/Raphael/Research/scripts/water_wire_viewing.tcl


  mol representation Tube 0.100000 6.000000
  mol color Name
  mol selection {protein}
  mol material Transparent
  mol addrep top



set input "results.txt"

set fid [open $input r]
set end_of_file 1
set line_counter 0
while {$end_of_file} {
  mol delrep 1 top 
  incr line_counter
  animate goto $line_counter
  set line [gets $fid]
  #puts "reading line $line_counter with [llength $line] elements"
  # this addes the helix ribbons




  # now assign the b value
  mol representation DynamicBonds 3.10000 0.300000 6.000000
  mol color Name
  mol selection index $line
  mol material Opaque
  mol addrep top
  mol selupdate 1 top 0
  after 200
  display update  
  if {[eof $fid]} {
    close $fid
    set end_of_file 0
    break
  }
}
