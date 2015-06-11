proc rotate_molecule {axis angle step} {
	#Rotate a molecule by an angle and take a picture every step
	#USAGE : rotate_molecule <axis> <angle> <step>
	set nstep [expr $angle / $step]
	display update on
	for { set i 1 } { $i <= $nstep } { incr i } {
		rotate $axis by $step
		display update
		puts [format "take picture #%04d at %8.3f deg." $i [expr $i*$step]]
		take_picture
}

