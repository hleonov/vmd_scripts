for {set i 1} {$i<=14} {incr i 1} {
	#mol load pdb "c$i\_rep_0.1.pdb"
	mol load pdb "./center_$i\.0.pdb"
	mol addrep top
	mol modselect 0 top "protein or resname ACE NAC"
	mol modselect 1 top "resname AMA"
	mol modstyle 0 top NewRibbons
	mol modstyle 1 top Lines
	mol modcolor 1 top colorID 1
	mol modcolor 0 top colorID 20
}
