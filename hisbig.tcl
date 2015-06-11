source ~/data/scripts/Tcl/bigtrj.tcl
source ~/data/scripts/Tcl/waitfor.tcl
source ~/data/scripts/Tcl/com.tcl

proc end {} {
        puts "Exiting..."
        exit
}

proc bin {z} {
	global slw 
	set sli [expr int( floor($slw * $z) ) + 1 ]
	return $sli
}

proc his {frame} {

	global fid countO countH ref all com rad m0list zoff wfreq
  	global fid countO countH nfra slw zorga nsl
	global sel1 sel2 sel3 sel4 org1 org2 org3 org4 y231 y232 y233 y234 y1491 y1492 y1493 y1494 y23a y149a

	incr nfra
	puts "$frame $nfra"
   
   	$all frame $frame
	$all move [measure fit $com $ref weight mass]
   
	set cntlist {}
	lappend cntlist [measure center $sel1]
	lappend cntlist [measure center $sel2]
	lappend cntlist [measure center $sel3]
	lappend cntlist [measure center $sel4]

 	set orglist {}
	lappend orglist [measure center $org1]
	lappend orglist [measure center $org2]
	lappend orglist [measure center $org3]
	lappend orglist [measure center $org4]

	set tyrlist23 {}
	lappend tyrlist23 [center_of_mass $y231]
	lappend tyrlist23 [center_of_mass $y232]
	lappend tyrlist23 [center_of_mass $y233]
	lappend tyrlist23 [center_of_mass $y234]

	set tyrlist149 {}
	lappend tyrlist149 [center_of_mass $y1491]
	lappend tyrlist149 [center_of_mass $y1492]
	lappend tyrlist149 [center_of_mass $y1493]
	lappend tyrlist149 [center_of_mass $y1494]

 	foreach monomer $m0list {
 
 		set cnt  [lindex $cntlist $monomer ]
 		set zorg [lindex [lindex $orglist $monomer] 2 ]
		set y23  [lindex [lindex $tyrlist23 $monomer] 2]
		set y149 [lindex [lindex $tyrlist149 $monomer] 2]
 
 		set wat [atomselect top "water and same residue as ( (x - [lindex $cnt 0])^2 + (y - [lindex $cnt 1] )^2 < $rad^2 )" frame $frame]
 		set zlist [$wat get {z}]
		
		#puts [measure minmax $wat]
 
 		for {set i 0} {$i < [$wat num]/3} {incr i} {

 			set idO   [expr 3*$i]
 			set idH1  [expr $idO + 1]
 			set idH2  [expr $idO + 2]

 			set z     [expr $zoff + [lindex $zlist  $idO] ]
 			#set sli   [expr int($slw * round($z)) ]
 			set sli   [bin $z]
 			incr countO($sli) 

 			set z     [expr $zoff + [lindex $zlist $idH1] ]
 			#set sli   [expr int($slw * round($z))]
 			set sli   [bin $z]
 			incr countH($sli) 

 			set z     [expr $zoff + [lindex $zlist $idH2] ]
 			#set sli   [expr int($slw * round($z))]
 			set sli   [bin $z]
 			incr countH($sli) 

 		}
 		$wat delete
 	}
 	set  zorga [expr $zorga + $zorg]
 	set  y23a  [expr $y23a  + $y23] 
 	set  y149a [expr $y149a + $y149] 

	if { $frame % $wfreq == 0 }  {

		set zorga [expr $zorga/double($nfra)]
		set y23a  [expr $y23a/double($nfra)]
		set y149a [expr $y149a/double($nfra)]

		for {set i -$nsl} {$i < $nsl} {incr i} {
			set $countO($i) [ expr double( $countO($i) )/double($nfra) ]
			set $countH($i) [ expr double( $countH($i) )/double($nfra) ]
			puts $fid "[expr double($i)/$slw] $countO($i) $countH($i) $zorga $y23a $y149a"
			flush $fid
		}  
	
		puts $fid " "
	
		puts "wrote; zeroring arrays and accumulators"

		set zorga 0.0
 		set y23a  0.0
 		set y149a 0.0
		set nfra  0

		for {set i -$nsl} {$i < $nsl} {incr i} {
			set countO($i) 0
			set countH($i) 0
		}
	}

}

set nsl 1000
set slw [expr 1/.1]
for {set i -$nsl} {$i < $nsl} {incr i} {
	set countO($i) 0
	set countH($i) 0
}

set wfreq [lindex $argv 0]
set rad   6.0
set zoff  0.0

set nfra  0
set zorga 0.0
set y23a  0.0
set y149a 0.0
set fid [open "waterdistrib_write_${wfreq}.dat" w+]

set mlist  [list 1 2 3 4]
set m0list [list 0 1 2 3]

# load
set idir /data/desrad-nb-o/jensenm/Aqp0/2B6Orun/desmondNPT_p40614/NPT_Berendsen_early_equil_2fs/
set my_system ${idir}/pdb/Aqp0_POPE_frame_1_cleaned_w_segnames.pdb
mol load pdb $my_system

# selections
set all [atomselect top "all"]
set ref [atomselect top "name CA" frame 0]
set com [atomselect top "name CA"]

$all moveby [ vecinvert [measure center $all] ]

puts "Original center: [measure center $all]"

foreach monomer $mlist {
	set sel${monomer}  [atomselect top "segname PRO${monomer} and resid 187 172 48" ]
	set bnd${monomer}  [atomselect top "segname PRO${monomer} and ( (resid 187 and name NE) or (resid 66 and name O) )"]
	set org${monomer}  [atomselect top "segname PRO${monomer} and ( (resid 68 184 and name ND2) )"]
	set y23${monomer}  [atomselect top "segname PRO${monomer} and ( (resid 23 and sidechain) )"]
	set y149${monomer} [atomselect top "segname PRO${monomer} and ( (resid 149 and sidechain) )"]
}

set lastfile 56172-65142.xtc
waitfor $lastfile end

cd /data/desrad-nb-o/jensenm/Aqp0/2B6Orun/rerun/xtc

#set lastfile foo
#bigdcd xtc his foo.xtc

bigdcd xtc his 0-11868.xtc 11868-45042.xtc 45048-56166.xtc 56172-65142.xtc

