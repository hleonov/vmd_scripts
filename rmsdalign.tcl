# Aligns whole system with molecule id $mol.
# That is selection $all (full system) is moved
# according to the transformation matrix computed
# from alignment of the selection $sel relative 
# to the reference coordinates of $sel at t=0 
#
# Example:
# Align system based on CA coordinates 
#
# "compile"
# source align.tcl
# Alternatively: source the script during vmd startup by
# introducing the line above into the file vmd.rc 
# placed in your home dir.. Hence, 
# the proc will always be available
#
proc align {sel ref all mol} {
   for {set frame 0} {$frame < [molinfo $mol get numframes]} {incr frame} {
      $sel frame $frame
      $all frame $frame
      $all move [measure fit $sel $ref]
   }
}

# Create selections
set all [atomselect top "all"]
set sel [atomselect top "protein and type CA"]
set ref [atomselect top "protein and type CA" frame 0]
#
# Execute 
#
align $sel $ref $all top


