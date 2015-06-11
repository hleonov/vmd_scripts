proc read_gro_xtc {gro xtc dt } {
	mol load gro $gro
	animate read xtc $xtc waitfor all
}

# define 3 reps: protein , dmpc, water (off)
proc def_display { } {
	mol addrep top
	mol addrep top
	mol modselect 0 top "protein or resname ACE NAC"
	mol modselect 1 top "resname DMPC"
	mol modselect 2 top "resname SOL"
	mol modstyle 0 top NewCartoon
	mol modstyle 1 top Lines
	mol modstyle 2 top Lines
	mol showrep top 2 off
	mol modcolor 1 top Name
	mol modcolor 2 top Name
}
