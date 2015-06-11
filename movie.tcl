proc take_picture {args} {
	global take_picture

	# when called with no parameter, render the image
	if {$args == {}} {

		#Take a picture every modulo frame
		if { [expr $take_picture(frame) % $take_picture(modulo)] == 0 }	{
		
			#IMAGE RENDERING ======================================= 
			#render $take_picture(method) $take_picture(tmpfile) $take_picture(options)\
			#       <$take_picture(tmpfile)\
			#       >[format "$take_picture(template)" $take_picture(frame)]
      			render $take_picture(method) $take_picture(tmp).dat \
			       [format "%s; convert -quality 100 -antialias %s.tga $take_picture(template)" \
			       $take_picture(options) $take_picture(tmp) $take_picture(frame)] 
		}

		# increase the count by one
		incr take_picture(frame)
		return
	}

	 lassign $args arg1 arg2
	# reset the options to their initial stat
	# (remember to delete the files yourself

	if {$arg1 == "reset"} {
		set take_picture(h) 240
		set take_picture(w) 320
		set take_picture(frame)  0
		set take_picture(tmp) "/Users/hleonov/tmp/vmd/tmp.r3d"
		set take_picture(template) "img"
		set take_picture(method) Raster3D
		set take_picture(options) "render -size $take_picture(w)x$take_picture(h) -aa -jpeg -quality 90"
		set take_picture(modulo) 1
		set take_picture(bindir) "/bioinfo/pbiladm/pbil/bin/vmd1.8.3/lib/vmd"
		return
	}


	# set one of the parameters
	if [info exists take_picture($arg1)] {
		if {[llength $args] == 1} {return "$arg1 is $take_picture($arg1)"}
		set take_picture($arg1) $arg2
		return
	}

	# otherwise, there was an error
	error {take_picture: [|reset|frame|format|method|modulo|h|w|option|tmp|template|bindir]}
	return
}

# to complete the initialization, this must be the first function
# called.  Do so automatically.
take_picture reset

proc make_trajectory_movie {} {

	set nframe [molinfo top get numframes]
	for {set i 0} {$i < $nframe} {incr i} {
		# go to the given frame
		animate goto $i
                #CENTER THE SYSTEM ON ITS CENTER OF MASS ===============
                #set sel [atomselect top "all"]
                #set center [list [measure center $sel]]
                #foreach mol [molinfo list] {
                #	molinfo $mol set center $center
                #	translate to [expr 0+$x] [expr 0+$y] [expr 0+$z]
                #}

                # force display update
                display update 
		# take the picture
		take_picture
	}
	return
}
