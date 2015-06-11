# this script will take as an inpust a selection (as a string)
# and then will move the entire molecues such that the center
# of the selection is at 0 0 0
#
# example:
# source ~/data/vmd_scripts/center_by_selection.tcl (not needed since it is the vmd.rc file)
# center_by_selection "resid 155 156 and protein"

proc center_by_selection {sel} {
  set selection [atomselect top $sel]
  set all [atomselect top "all"]
  puts [measure center $selection]
  $all moveby [vecinvert [measure center $selection]]
}
