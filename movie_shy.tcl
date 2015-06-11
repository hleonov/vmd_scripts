
#This will generate a movie out of a trajectory 
#skipping frames as specified in freq
#rotating the molecule by "angle" degrees during the whole movie, around "axis".
#output: *.pov
proc make_movie_rot_trj {freq angle axis} {
   # get the number of frames in the movie
   set num [molinfo top get numframes]
   set step [expr $angle/[expr $num/"$freq.0"]]
   puts $step
   # loop through the frames
   for {set i 0} {$i < $num} {incr i $freq} {
      # go to the given frame
      animate goto $i
      # for the display to update
      display update
      # take the picture
	  rotate $axis by $step
      set filename name.[format "%04d" [expr $i/$freq]].pov
      render POV3 $filename
   }
}

#This will generate a movie out of a trajectory 
#skipping frames as specified in freq
#output: *.pov
proc make_trajectory_movie {freq} {
   # get the number of frames in the movie
   set num [molinfo top get numframes]
   # loop through the frames
   for {set i 0} {$i < $num} {incr i $freq} {
      # go to the given frame
      animate goto $i
      # for the display to update
      display update
      # take the picture
      set filename name.[format "%04d" [expr $i/$freq]].pov
      render POV3 $filename
   }
}
