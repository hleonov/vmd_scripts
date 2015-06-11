set fid [open "count_waters.dat" w]
set waters [atomselect top "resname HOH SOL and type OW"];
foreach ow [$waters list] {
	set ow_sel [atomselect top "index $ow"]
	set rid [$ow_sel get resid]
	#puts "atom: $ow\t residue: $rid"
	set sel [atomselect top "type OW and (within 3.5 of index $ow)"]
	puts $fid [$sel num] 
}

close $fid



