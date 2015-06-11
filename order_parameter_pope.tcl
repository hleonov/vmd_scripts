source bigdcd.tcl

proc order { frame } {

    global fid1 fid2 fid3 fid4

    puts $frame

    # c18 monolayer 1 AKA leg 1

    for { set i 2 } { $i < 19 } { incr i } {
    	set od1($i) 0.0
    	if {$i == 9 || $i == 10} { 
    	    set c [ atomselect top "resname POPE and (name C2${i} H${i}1 )" frame $frame ]
    	} elseif { $i == 18 } {
    	    set c [ atomselect top "resname POPE and (name C2${i} H${i}R H${i}S H${i}T)" frame $frame ]
    	} else { 
    	    set c [ atomselect top "resname POPE and (name C2${i} H${i}R H${i}S )" frame $frame ]
    	}
    	set idc [$c list]
    	set crd [$c get {x y z}]
    	$c delete
    	
    	# all cases but double bond and terminus
    	
    	if { $i != 9 || $i != 10 || $i != 18 } {	    
    	    for { set j 0 } { $j < [expr [llength $idc] / 3 ] } { incr j } {
    		set cc      [ expr $j*3 ]
    		set cch1    [ expr ( $j*3 ) + 1 ]
    		set cch2    [ expr ( $j*3 ) + 2 ]
    		set crdC    [ lindex $crd $cc ]
    		set crdH1   [ lindex $crd $cch1 ]
    		set crdH2   [ lindex $crd $cch2 ]
    		set r1      [ vecsub $crdH1 $crdC ]
    		set r2      [ vecsub $crdH2 $crdC ]	
    		set r1l     [ veclength $r1 ]
    		set r2l     [ veclength $r2 ]
    		set r1l     [ expr 1.0 / $r1l ]
    		set r2l     [ expr 1.0 / $r2l ]
    		set r1      [ vecscale $r1l $r1 ]
    		set r2      [ vecscale $r2l $r2 ]
    		set o1      [ lindex $r1 2]
    		set o1      [ expr 1.5 * ($o1 * $o1) - 0.5 ]
    		set o2      [ lindex $r2 2]
    		set o2      [ expr 1.5 * ($o2 * $o2) - 0.5 ]
    		set od1($i) [ expr $od1($i) + ( 0.5 * ($o1 + $o2) ) ] 
    	    }
    	}
    	
    	# case double bond
    	
    	if { $i == 9 || $i == 10 } {
    	    for { set j 0 } { $j < [expr [llength $idc] / 2 ] } { incr j } {
    		set cc      [ expr $j*2 ]
    		set cch1    [ expr ( $j*2 ) + 1 ]
    		set crdC    [ lindex $crd $cc ]
    		set crdH1   [ lindex $crd $cch1 ]
    		set r1      [ vecsub $crdH1 $crdC ]
    		set r1l     [ veclength $r1 ]
    		set r1l     [ expr 1.0 / $r1l ]
    		set r1      [ vecscale $r1l $r1 ]
    		set o1      [ lindex $r1 2]
    		set o1      [ expr 1.5 * ($o1 * $o1) - 0.5 ]
    		set od1($i) [ expr $od1($i) + $o1 ] 
    	    }
    	}
    	
    	# case terminus
    	
    	if { $i == 18 } {
    	    for { set j 0 } { $j < [expr [llength $idc] / 4 ] } { incr j } {
    		set cc      [ expr $j*3 ]
    		set cch1    [ expr ( $j*3 ) + 1 ]
    		set cch2    [ expr ( $j*3 ) + 2 ]
    		set cch3    [ expr ( $j*3 ) + 3 ]
    		set crdC    [ lindex $crd $cc ]
    		set crdH1   [ lindex $crd $cch1 ]
    		set crdH2   [ lindex $crd $cch2 ]
    		set crdH3   [ lindex $crd $cch3 ]
    		set r1      [ vecsub $crdH1 $crdC ]
    		set r2      [ vecsub $crdH2 $crdC ]
    		set r3      [ vecsub $crdH3 $crdC ]
    		set r1l     [ veclength $r1 ]
    		set r2l     [ veclength $r2 ]
    		set r3l     [ veclength $r3 ]
    		set r1l     [ expr 1.0 / $r1l ]
    		set r2l     [ expr 1.0 / $r2l ]
    		set r3l     [ expr 1.0 / $r3l ]
    		set r1      [ vecscale $r1l $r1 ]
    		set r2      [ vecscale $r2l $r2 ]
    		set r3      [ vecscale $r3l $r3 ]
    		set o1      [ lindex $r1 2]
    		set o1      [ expr 1.5 * ($o1 * $o1) - 0.5 ]
    		set o2      [ lindex $r2 2]
    		set o2      [ expr 1.5 * ($o2 * $o2) - 0.5 ]
    		set o3      [ lindex $r3 2 ]
    		set o3      [ expr 1.5 * ($o3 * $o3) - 0.5 ]
    		set od1($i) [ expr $od1($i) + ( 0.33 * ( $o1 + $o2 + $o3) ) ]
    	    }
    	}
    }

    # write average for each frame

    for { set id 2 } { $id  < 19 } { incr id } {
    	if {$id == 9 || $id == 10} {
    	    puts $fid1($is) "[expr $od1($id) / ( [llength $idc] / 2 ) ]"
    	    flush $fid1($is)
    	} elseif {$id == 18} {
    	puts $fid1($is) "[expr $od1($id) / ( [llength $idc] / 4 ) ]"
    	    flush $fid1($is)
    	} else {
    	    puts $fid1($is) "[expr $od1($id) / ( [llength $idc] / 3 ) ]"
    	    flush $fid1($is)
    	}
    }  

    # c16 monolayer 1 AKA leg 2

    for { set i 2 } { $i < 17 } { incr i } {
    	set od1($i) 0.0
    	if { $i == 16 } {
    	    set c [ atomselect top "resname POPE and (name C2${i} H${i}X H${i}Y H${i}Z) and 
    		    (same residue as name C31 and z > 0 and 
    		    (not within $rmin of protein or within $rmin of resname ETAM) and 
    		    (within $rmax of protein or within $rmax of resname ETAM))" frame $frame ]
    	} else { 
    	    set c [ atomselect top "resname POPE and (name C2${i} H${i}Y H${i}Y ) and 
    		    (same residue as name C31 and z > 0 and 
    		    (not within $rmin of protein or within $rmin of resname ETAM) and 
    		    (within $rmax of protein or within $rmax of resname ETAM))" frame $frame ]
    	}
    	set idc [$c list]
    	set crd [$c get {x y z}]
    	$c delete   
    	
    	# case not terminus 
    	
    	if { $i != 16 } {
    	    for { set j 0 } { $j < [expr [llength $idc] / 3 ] } { incr j } {
    		set cc      [ expr $j*3 ]
    		set cch1    [ expr ( $j*3 ) + 1 ]
    		set cch2    [ expr ( $j*3 ) + 2 ]
    		set crdC    [ lindex $crd $cc ]
    		set crdH1   [ lindex $crd $cch1 ]
    		set crdH2   [ lindex $crd $cch2 ]
    		set r1      [ vecsub $crdH1 $crdC ]
    		set r2      [ vecsub $crdH2 $crdC ]	
    		set r1l     [ veclength $r1 ]
    		set r2l     [ veclength $r2 ]
    		set r1l     [ expr 1.0 / $r1l ]
    		set r2l     [ expr 1.0 / $r2l ]
    		set r1      [ vecscale $r1l $r1 ]
    		set r2      [ vecscale $r2l $r2 ]
    		set o1      [ lindex $r1 2]
    		set o1      [ expr 1.5 * ($o1 * $o1) - 0.5 ]
    		set o2      [ lindex $r2 2]
    		set o2      [ expr 1.5 * ($o2 * $o2) - 0.5 ]
    		set od1($i) [ expr $od1($i) + ( 0.5 * ($o1 + $o2) ) ] 
    	    }
    	}
    	
    	# case terminus
    	
    	if { $i == 16 } {
    	    for { set j 0 } { $j < [expr [llength $idc] / 4 ] } { incr j } {
    		set cc      [ expr $j*3 ]
    		set cch1    [ expr ( $j*3 ) + 1 ]
    		set cch2    [ expr ( $j*3 ) + 2 ]
    		set cch3    [ expr ( $j*3 ) + 3 ]
    		set crdC    [ lindex $crd $cc ]
    		set crdH1   [ lindex $crd $cch1 ]
    		set crdH2   [ lindex $crd $cch2 ]
    		set crdH3   [ lindex $crd $cch3]
    		set r1      [ vecsub $crdH1 $crdC ]
    		set r2      [ vecsub $crdH2 $crdC ]
    		set r3      [ vecsub $crdH3 $crdC ]
    		set r1l     [ veclength $r1 ]
    		set r2l     [ veclength $r2 ]
    		set r3l     [ veclength $r3 ]
    		set r1l     [ expr 1.0 / $r1l ]
    		set r2l     [ expr 1.0 / $r2l ]
    		set r3l     [ expr 1.0 / $r3l ]
    		set r1      [ vecscale $r1l $r1 ]
    		set r2      [ vecscale $r2l $r2 ]
    		set r3      [ vecscale $r3l $r3 ]
    		set o1      [ lindex $r1 2]
    		set o1      [ expr 1.5 * ($o1 * $o1) - 0.5 ]
    		set o2      [ lindex $r2 2]
    		set o2      [ expr 1.5 * ($o2 * $o2) - 0.5 ]
    		set o3      [ lindex $r3 2 ]			    
    		set o3      [ expr 1.5 * ($o3 * $o3) - 0.5 ]
    		set od1($i) [ expr $od1($i) + ( 0.33 * ( $o1 + $o2 + $o3) ) ]
    	    }
    	}
    }

    # write data for monolayer 1 leg 2 to file

    for { set id 2 } { $id  < 17 } { incr id } {
    	if {$id != 17} {
    	    puts $fid3($is) "[expr $od1($id) / ( [llength $idc] / 3 ) ]"
    	    flush $fid3($is)
    	} else {
    	    puts $fid3($is) "[expr $od1($id) / ( [llength $idc] / 4 ) ]"
    	    flush $fid3($is)
    	}
    }
}

set sys POPE

    set name1 "ord_${sys}_rad_c18_1_${is}.out"
    set fid1($is) [open $name1 w+]    
    set name2 "ord_${sys}_rad_c18_2_${is}.out"
    set fid2($is) [open $name2 w+]    
 
set inp gA_POPE

mol load psf ${sys}/${inp}.psf
animate read pdb ${sys}/${inp}.pdb

cd dcd

bigdcd order POPE/output2/gA_POPE.md.2.dcd POPE/output3/gA_POPE.md.3.dcd POPE/output4/gA_POPE.md.4.dcd POPE/output5/gA_POPE.md.5.dcd POPE/output6/gA_POPE.md.6.dcd POPE/output7/gA_POPE.md.7.dcd

