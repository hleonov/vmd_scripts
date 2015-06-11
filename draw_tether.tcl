# this script will operate on the top molecule
# it assumed a trajectory is loaded already
# it will also align all frame to the original frame according tto the protien Ca
#
# It will then draw the position of the tether point as a sphere.
# The position of the tether point is specified by the time and velocity when computing z
# and the x and y are taken from the original poistion and the x and y of the selection ceinter
# thus we draw 2 points.
#
# the inputs are:
# 1. The time difference between each frame ($dt)
# 2. The velocity of the tether point ($v)
# 3. The selection of which you will take the x and y from ($sel)

proc draw_tether {dt v sel} {
  # selections for aligment
  set ref [atomselect top "protein and name CA" frame 0]
  set all [atomselect top "all"]
  set ca  [atomselect top "protein and name CA"]
  # get the x, y and z of the original position of the target
  set x0 [lindex [measure center $sel weight mass] 0]
  set y0 [lindex [measure center $sel weight mass] 1]
  set z0 [lindex [measure center $sel weight mass] 2]
    graphics top delete all
    # now run throu the frames
    for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame} {
      # update the display
      animate goto $frame
      display update
      # update selections
      $ca frame $frame
      $all frame $frame
      $sel frame $frame
      # align it to the original
      $all move [measure fit $ca $ref]
      # get the x and y of the current position of the target
      set x [lindex [measure center $sel weight mass] 0]
      set y [lindex [measure center $sel weight mass] 1]
      # calculate the z positions
      set z [expr ($v * $dt * $frame) + $z0]
      # remove the following line if you want to have contiuous drawing
      graphics top delete all
      graphics top material Transparent
      graphics top color 2
      draw sphere  [list $x $y $z ] radius 0.5
      graphics top color 8
      draw sphere  [list $x0 $y0 $z ] radius 0.5
   }
}





