
proc view_cylinder {sel radius length} {
  set center [measure center $sel]
  draw cylinder [list [lindex $center 0] [lindex $center 1] [expr [lindex $center 2]-$length]] [list [lindex $center 0] [lindex $center 1] [expr [lindex $center 2]+$length]] radius $radius

}

#graphics top delete all

proc draw_it {sel radius clr length} {
  graphics top color $clr
  graphics top material Transparent
  view_cylinder $sel $radius $length
}

