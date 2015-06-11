set nac [atomselect top "resname NAC"]
set ace [atomselect top "resname ACE"]
set center_nac [measure center $nac]
set center_ace [measure center $ace]
set prot_vec [vecsub $center_ace $center_nac]
set cos_theta [vecdot [vecnorm $prot_vec] {0 0 1}]
set theta [expr {acos($cos_theta)}]
set theta_deg [expr {$theta*180/3.14159265359}]



proc angle_between_vectors_in_radians {vec1 vec2} {
        set cos_of_angle [vecdot [vecnorm $vec1] [vecnorm $vec2]]

        # if cos of angle is a little bit more
        # than 1, round it to 1 to avoid an error.
        # If it is significantly more than 1, complain.
        if {$cos_of_angle>1} {
                if {$cos_of_angle<1.0000001} {
                        set cos_of_angle 1.0
                } else {
                        puts "ERROR: cannot take arccos of $cos_of_angle"
                                # tcl also outputs an error message, but
                                # this one makes it clear what went wrong
                }
        }
        if {$cos_of_angle<-1} {
                if {$cos_of_angle>-1.0000001} {
                        set cos_of_angle -1.0
                } else {
                        puts "ERROR: cannot take arccos of $cos_of_angle"
                }
        }

        return [expr {acos($cos_of_angle)}]

}

proc radians_to_degrees {angle_in_radians} {
        # pi = acos(-1) = 3.14159265359
        return [expr {$angle_in_radians*180/3.14159265359}]
} 
