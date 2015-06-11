#Fix bonds in POP(E) bilayer in order to proper displaying in VMD

proc fixbonds {mol} {

	
	set O7  [atomselect $mol "resname \"POP.*\" and name  O7"]
	set P8  [atomselect $mol "resname \"POP.*\" and name  P8"]
	set O11 [atomselect $mol "resname \"POP.*\" and name O11"]
	
	foreach O [$O7 list] P [$P8 list] OO [$O11 list] {
		#07
		set sel [atomselect $mol "index $O"]
		set blist {}
		lappend blist [expr $O - 1] [expr $O + 1]
		$sel setbonds [list $blist]
		$sel delete
		#P8
		set sel [atomselect $mol "index $P"]
		set blist {}
		lappend blist  [expr $P - 1] [expr $P + 1] [expr $P + 2] [expr $P + 3]
		$sel setbonds [list $blist]
		$sel delete
		#011
		set sel [atomselect $mol "index $OO"]
		set blist {}
		lappend blist [expr $OO - 3] [expr $OO + 1]
		$sel setbonds [list $blist]
		$sel delete
	}
}
#fixbonds top
