package require Orient
#namespace import Orient::orient

proc compute_pc {seltxt} {
	set sel [atomselect top "$seltxt"]
	set I [draw principalaxes $sel]
	namespace forget Orient::orient
	return $I
}

proc norm_to_plane {v1 v2} {
	set cross [vecnorm [veccross $v2 $v1]]
	return $cross
}


