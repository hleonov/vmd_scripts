
proc view_cylinder {sel radius} {
  set center [measure center $sel]
  draw cylinder [list [lindex $center 0] [lindex $center 1] [expr [lindex $center 2]-20]] [list [lindex $center 0] [lindex $center 1] [expr [lindex $center 2]+20]] radius $radius

}
#set sf1 [atomselect top "resid 155 and protein"]

graphics top delete all

proc draw_it {sel radius step} {
  set cintv [expr 1040 - 17]
  set n [molinfo top get numframes]
  for {set i 0} {$i < $n} {incr i $step} {
    set clr [expr 17+ $i*$cintv/$n]
    graphics top color $clr
    $sel frame $i
    view_cylinder $sel $radius
    puts $clr

  }
}



