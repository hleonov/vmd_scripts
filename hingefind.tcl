####--------------------------------------------------------------------
####------  H I N G E F I N D  - TCL PLUG-IN FOR VMD 5/30/97 -----------
####--------------------------------------------------------------------
#   this software is copyrighted, (c) 1995-1997, by Willy Wriggers 
#   under the terms of the legal statement in the distribution 
#   (ftp://ftp.scripps.edu/pub/wriggers/hingefind/legal).
#
#### DOCUMENTATION -----------------------------------------------------
#    web site: ftp://ftp.scripps.edu/pub/wriggers/hingefind/hingefind.html
#
#### VMD ---------------------------------------------------------------
#    this program requires the free molecular graphics package "vmd".
#    web site: http://www.ks.uiuc.edu/Research/vmd
#
#### LIST OF GLOBAL VARIABLES ------------------------------------------
#    all1             list of all atom indices (molecule 1)
#    ca1              list of search atom indices (molecule 1)
#    ca2              list of search atom indices (molecule 2)
#    cutdom           minimum number of atoms per domain
#    dindex           list of domain membership indices
#    domain_count     number of found domains
#    domindex1        array of domain atom indices (molecule 1)
#    domindex2        array of domain atom indices (molecule 2)
#    full1            list of comparison atom indices (molecule 1)
#    full2            list of comparison atom indices (molecule 2)
#    hinge_ref_id     id of current reference domain for hinge-drawing
#    initrad          initial radius of seed subset in A
#    loadreps1        number of representations of molecule 1
#    maxcounter       maximum number of iterations per domain
#    mol1             id of molecule 1
#    mol2             id of molecule 2
#    ndomains         maximum possible number of domains
#    sup_rmsd         list of best rmsds (rmsd history) 


#### LOAD -------------------------------------------------------------- 
proc load {} {
#   - load the structures, select atom subsets and set parameters.
#   - this procedure should be edited by the user and is called once 
#     at the beginning of a session.

    global all1 ca1 ca2 mol1 mol2
    global full1 full2 
    global sup_rmsd dindex
    global initrad maxcounter ndomains cutdom
    global hinge_ref_id loadreps1

    # initial radius of seed subset in A; should be > 15.
    set initrad 20.0
    # maximum number of iterations until a domain has converged.
    set maxcounter 20
    # maximum possible number of domains; set large for full search.
    set ndomains 999
    # minimum number of atoms which a domain should have (> 4).
    set cutdom 10

    # load the first pdb file from pdb databank (pdb entry 1lfg)
    mol load gro from_EM_bfgs.gro
    # (or read a file with "mol load pdb foo.pdb"!)
    set mol1 [molinfo top]

    # load the second pdb file from pdb databank (pdb entry 1lfh)
    mol load pdb model1_ref.pdb
    # (or read a file with "mol load pdb bar.pdb"!)
    set mol2 [molinfo top]

    # customize the rendering of the molecules here:
    mol delrep 0 $mol1
    mol representation Tube 0 1
    mol color ColorID 11 
    mol selection {all}
    mol addrep $mol1
    mol delrep 0 $mol2
    mol representation Tube 0 1
    mol color ColorID 8
    mol selection {all}
    mol addrep $mol2


    set all1 [atomselect $mol1 all]
    $all1 global

    # select the compared atoms in the first molecule 
    set full1 [atomselect $mol1 "(resid 1 to 295) and type CA"]
    $full1 global
    set ca1 $full1

    # select the compared atoms in the second molecule 
    set full2 [atomselect $mol2 "type CA"]
    $full2 global
    set ca2 $full2

    # check input
    if {[$ca1 num] != [$ca2 num]} {
	error "load> two selections need same number of atoms!"
    }   
    if {$cutdom < 5} {
	error "load> cutdom needs to be larger than 4!"
    }   

    # the following is necessary to determine the rendering state
    # do not change
    # indicate that no hinge axes are drawn
    set hinge_ref_id -1
    # remember number of representations of molecule 1
    set loadreps1 [molinfo $mol1 get numreps]

    return
}


#### VECDIST ----------------------------------------------------------- 
proc vecdist {v1 v2} {
#   - general purpose fast 3dim vector distance procedure.
#   - factor 2 speedup relative to vecnorm [vecsub ...]

    lassign $v1 x1 x2 x3
    lassign $v2 y1 y2 y3
    return [expr sqrt(pow($x1-$y1,2) + pow($x2-$y2,2) + pow($x3-$y3,2))] 
}


#### LSUBTRACT --------------------------------------------------------- 
proc lsubtract {l1 l2} {
#   - general purpose list subtraction.
#   - find those atoms in l1 which are not in l2.

    set result {}
    foreach element $l1 {
	if {[lsearch $l2 $element] == -1} {
	    lappend result $element
	}
    }
    return $result
}


#### CRITERION --------------------------------------------------------- 
proc criterion {sel1 sel2 eps} {
#   - find those atoms which have rmsd < eps in the two selections
#   - return the info as a list of lists of indicies
#   - called by procedure 'superimpose'.

    set coord1 [$sel1 get {x y z}]
    set coord2 [$sel2 get {x y z}]
    set atomid1 [$sel1 get index]
    set atomid2 [$sel2 get index]
    set result1 {}
    set result2 {}

    foreach c1 $coord1 c2 $coord2 id1 $atomid1 id2 $atomid2 {
	if { [vecdist $c2 $c1] < $eps} {
	    lappend result1 $id1
	    lappend result2 $id2
	}
    }
    set result [list $result1 $result2] 
    return $result
}


#### SUPERIMPOSE ------------------------------------------------------- 
proc superimpose {eps} {
#   - superposition based on the atoms which match the best.
#   - return the best matching subsets as a list of lists of indicies.
#   - called by procedure 'convergence'.
    
    global all1 ca1 ca2 mol1 mol2

    # get atoms that are within 'eps'
    set survivor [criterion $ca1 $ca2 $eps]
    set survivor1 [lindex $survivor 0]
    set survivor2 [lindex $survivor 1]
    if {[llength $survivor1] < 5} {
       return $survivor
    }

    # and fit by the best matching set
    set survivor_sel1 [atomselect $mol1 "index $survivor1"]
    set survivor_sel2 [atomselect $mol2 "index $survivor2"]
    set t [measure fit $survivor_sel1 $survivor_sel2 weight mass]
    $all1 move $t
    display update
    return $survivor
}


#### SEED -------------------------------------------------------------- 
proc seed {sel1 sel2 radius} {
#   - find a seed subset (all atoms of 'sel1' within 'radius' of the
#     first atom in 'sel1').
#   - return the info as a list of lists of indicies.
#   - called by procedure 'convergence'.

    set coord1 [$sel1 get {x y z}]
    set atomid1 [$sel1 get index]
    set atomid2 [$sel2 get index]

    set cd1start [lindex $coord1 0]
    set id1start [lindex $atomid1 0]

    set result1 {}
    set result2 {}
    foreach cd1 $coord1 id1 $atomid1 id2 $atomid2 {
	if { [vecdist $cd1 $cd1start] < $radius} {
	    lappend result1 $id1
	    lappend result2 $id2
	}
    }
    set result [list $result1 $result2] 
    return $result
}


#### CONVERGENCE ------------------------------------------------------- 
proc convergence {eps radius} {
#   - iterate until a domain is found.
#   - return domain as a list of lists of atom indices.  
#   - called by procedure 'partition'.

    global all1 ca1 ca2 mol1 mol2 maxcounter

    # initiate search in seed-subset
    set startids [seed $ca1 $ca2 $radius]
    set startids1 [lindex $startids 0]
    set startids2 [lindex $startids 1]
    if {[llength $startids1] < 5} {
       # nothing to find in seed set, return formally as a domain
       return $startids
    }

    # if seed-subset sufficiently large, proceed with initial fit
    set startsel1 [atomselect $mol1 "index $startids1"]
    set startsel2 [atomselect $mol2 "index $startids2"]
    set t [measure fit $startsel1 $startsel2 weight mass]
    $all1 move $t
    
    # search iteratively for the domain 
    set count 0
    set prev [superimpose $eps]
    set curr [superimpose $eps]
    while {[string compare $prev $curr]} {
	set prev $curr
	set curr [superimpose $eps]	
	incr count
	if {$count == $maxcounter} {
	    puts "convergence - warning: a domain did not converge."
	    break
	}
    }
    if {[llength [lindex $curr 0]] < 5} {
       # even though seed set was sufficiently large, 
       # the found domain is too small. return nothing and 
       # let 'partition' try again with smaller radius
       return {{} {}}
    }
    return $curr
}


#### UPDATE_BOUNDARIES ------------------------------------------------- 
proc update_boundaries {eps number} {
#   - update the 'sup_rmsd' and 'dindex' lists with the 'number' domain.
#   - actually "save" the domain and update it's boundaries with earlier
#     found domains.
#   - called by procedure 'partition'.

    global full1 full2 
    global sup_rmsd dindex

    set coord1 [$full1 get {x y z}]
    set coord2 [$full2 get {x y z}]

    set newsup_rmsd {}
    set newdindex {}

    foreach cd1 $coord1 cd2 $coord2 su $sup_rmsd di $dindex {      
        set dist [vecdist $cd1 $cd2] 
	if {$dist < $su} {
	    lappend newsup_rmsd $dist
            if {$dist < $eps} {
		lappend newdindex $number
            } else {
		lappend newdindex $di
            }
        } else {
            lappend newsup_rmsd $su
            lappend newdindex $di
        }
    }

    set sup_rmsd $newsup_rmsd
    set dindex $newdindex

    return
}


#### SORT_N_RENDER ----------------------------------------------------- 
proc sort_n_render {} {
#   - sort the found domains, color-code and render, output info, 
#     update 'domain_count'.
#   - called by procedure 'partition'.

    global domindex1 domindex2
    global full1 full2 all1
    global mol1 
    global dindex cutdom domain_count
    
    # bubblesort
    puts "sort_n_render> now sorting the domains by size."
    for {set i 1} {$i < $domain_count} {incr i} {
    for {set j [expr $domain_count - 1]} {$j >= $i} {incr j -1} {
	set j1 [expr $j - 1]
	if {[llength $domindex1($j1)] < [llength $domindex1($j)]} {
	set xchange $domindex1($j1)
	set domindex1($j1) $domindex1($j)
	set domindex1($j) $xchange
	set xchange $domindex2($j1)
	set domindex2($j1) $domindex2($j)
	set domindex2($j) $xchange
    } }	}

    # render sorted domains by color-coded residues
    set t [measure fit $full1 $full2 weight mass]
    $all1 move $t
    display update
    set ntot 0
    set nrtot 0
    for {set i 0} {$i < $domain_count} {incr i} {
        set nct [llength $domindex1($i)] 
	if {$nct < $cutdom} {
		break
	}
	set ntot [expr $ntot + $nct]

        mol rep Tube 0.3 6
        mol color colorid [expr $i % 17]
        set sel [atomselect $mol1 "index $domindex1($i)"]
	# color by residue
        set rlist [$sel get residue] 
        set resindex {}
        foreach slist $rlist {
            lappend resindex [lindex $slist 0]
        }
        set reslist [lsort -integer $resindex]
	set nrtot [expr $nrtot + [llength $reslist]]
        mol selection "residue $reslist"
        mol addrep $mol1

        #output resids
        set rlist [$sel get resid] 
        set resindex {}
        foreach slist $rlist {
            lappend resindex [lindex $slist 0]
        }
        set reslist [lsort -integer $resindex]
        puts "sort_n_render> domain number: $i - number of atoms: $nct - resids:"     
	puts "$reslist"
    }
    # from now on only consider domains larger than 'cutdom' atoms
    puts "sort_n_render> total number of domains larger than $cutdom atoms: $i ."
    set domain_count $i
    puts "sort_n_render> number of atoms involved: $ntot ."
    if {$ntot != $nrtot } {
	puts "sort_n_render> warning: there are more than one atoms than residues, rendering by residue may not be accurate!"   
    }
    set rest [expr [$full1 num] - $ntot]
    puts "sort_n_render> unconverged rest: $rest atoms."
    return 
}


#### PARTITION --------------------------------------------------------- 
proc partition {eps} {
#   - the main partitioning procedure.
#   - generates and renders the domains.
#   - stores the domains of molecules 1 and 2 in 'domindex' arrays. 
#   - called by the user.

    global all1 ca1 ca2 mol1 mol2
    global full1 full2 
    global sup_rmsd dindex
    global initrad ndomains cutdom
    global domain_count domindex1 domindex2
    global hinge_ref_id loadreps1

    # delete any drawn hinge axes
    if {$hinge_ref_id != -1} {
	set hinge_ref_id -1 
	draw delete all
    }
    # reset domain representations
    set rendernum [molinfo $mol1 get numreps]
    while {$rendernum > $loadreps1} {
	mol delrep [expr $rendernum -1] $mol1
	set rendernum [molinfo $mol1 get numreps]
    }


    set ca1 $full1
    set ca2 $full2

    # initialize
    set sup_rmsd {}
    set dindex {}
    for {set i 0} {$i < [$full1 num]} {incr i} { 
        lappend sup_rmsd 999
        lappend dindex -1
    }
    set domain_count 0
    set seedradius $initrad 

    # partition the protein
    while {$domain_count < $ndomains} {

	# find a domain        
 	set domain [convergence $eps $seedradius]
        set domain1 [lindex $domain 0]
        set domain2 [lindex $domain 1]

	if {[llength $domain1] == 0} {
	    # convergence found nothing, try smaller seedradius
	    puts "partition> trying 20% smaller seed radius..."
            set seedradius [expr $seedradius * 0.8]
            continue 
	}		
        
	# remove from ca1 and ca2 those atoms which have been found
	set ca1_atoms [$ca1 list]
	set disjoint1 [lsubtract $ca1_atoms $domain1]
	set ca2_atoms [$ca2 list]
	set disjoint2 [lsubtract $ca2_atoms $domain2]
	
	# and reset the definitions of ca1 and ca2
	if {[llength $disjoint1] > 0} {
	    set ca1 [atomselect $mol1 "index $disjoint1"]
	    $ca1 global
	    set ca2 [atomselect $mol2 "index $disjoint2"]
	    $ca2 global
	}

	if {[llength $domain1] > 4} {
	    # convergence found a domain after normal iteration
            # reset seed radius and update domain boundaries
            set seedradius $initrad
            update_boundaries $eps $domain_count
	}	

	incr domain_count
	if {$domain_count == 1} {
		puts "partition> now have 1 domain."
	} else {
		puts "partition> now have $domain_count domains."
	} 

	if {[llength $disjoint1] < $cutdom} {
	    puts "partition> protein now fully partitioned."
	    break
	}	 
    }

    # store the unsorted domains in the 'domindex' arrays
    set domindex1(0) {}
    set domindex2(0) {}
    for {set i 1} {$i < $domain_count} {incr i} {
	set domindex1($i) {}
	set domindex2($i) {}
    }   
    set atomid1 [$full1 get index]
    set atomid2 [$full2 get index]
    foreach di $dindex id1 $atomid1 id2 $atomid2 {
    if {$di >= 0} { 
         lappend domindex1($di) $id1
         lappend domindex2($di) $id2
    } }

    sort_n_render 
    return
}


#### HINGE ------------------------------------------------------------- 
proc hinge {domid1 domid2} {
#   - generate and render the effective rotation axis of the movement 
#     of 'domid2' (moving domain) relative to 'domid1' (reference
#     domain).
#   - called by the user.

    global all1 mol1 mol2
    global full1 full2 
    global domain_count domindex1 domindex2
    global hinge_ref_id

    # set radius of drawn effect. rot. axis
    set raxis 0.5 
    # set length of drawn effect. rot. axis
    set laxis 80 

    # check if domain ids are within range
    if {($domid1 < 0)||($domid1 >= $domain_count)||($domid2 < 0)||($domid2 >= $domain_count)} {
	error "hinge> domain id out of range!"
    } 

    # reset drawn hinge axes if reference domain changed
    if {$domid1 != $hinge_ref_id} {
	set hinge_ref_id $domid1 
	puts "hinge> new reference domain, reset drawn hinge axes."
	draw delete all
    }
	
    puts "hinge> effective rotation evaluated for domain: $domid2 ."

    set selid1_1 [atomselect $mol1 "index $domindex1($domid1)"]
    set selid1_2 [atomselect $mol2 "index $domindex2($domid1)"]
    set selid2_1 [atomselect $mol1 "index $domindex1($domid2)"]
    set selid2_2 [atomselect $mol2 "index $domindex2($domid2)"]

    # superimpose the two structures by reference domain (domid1)
    set t [measure fit $selid1_1 $selid1_2 weight mass]
    $all1 move $t
    puts "hinge> structures superimposed by domain $domid1 (reference domain)."

    # compute the COM of the moving domain in the two structures
    set com1 [measure center $selid2_1 weight mass]
    set com2 [measure center $selid2_2 weight mass]
    set cdis [vecdist $com1 $com2]
    puts "hinge> com-movement of domain $domid2: $cdis ."

    # calculate bisecting point and normalvector of bisecting plane
    set bi [vecscale [vecadd $com1 $com2] 0.5]
    set pl [vecnorm [vecsub $com2 $com1]]

    # calculate the best fit rmsd of the moving domain and reset
    set trans_matrix [measure fit $selid2_1 $selid2_2 weight mass]
    $all1 move $trans_matrix
    set rmsid [measure rmsd $selid2_1 $selid2_2 weight mass]
    set t [measure fit $selid1_1 $selid1_2 weight mass]
    $all1 move $t

    # compute rotation matrix for rotations about com2
    set t21 [transoffset [vecsub $com1 $com2]]
    set rot2 [transmult $trans_matrix $t21]
    
    # rotate com1 twice about com2
    set p1 [coordtrans $rot2 $com1]
    set p2 [coordtrans $rot2 $p1]

    # compute the rotation axis from the 3 points
    set rideal [vecnorm [veccross [vecsub $com1 $p2] [vecsub $com1 $p1]]]
	
    # project com2 onto the 3-point-plane
    set new [vecsub $com2 [vecscale $rideal [vecdot $rideal [vecsub $com2 $com1]]]]

    # compute rotation angle (least squares fit)
    set cosine [vecdot [vecnorm [vecsub $new $com1]] [vecnorm [vecsub $new $p1]]]
    set angl [expr acos($cosine)]

    # compute projection of rot axis on bisecting plane
    set perp [vecdot $rideal $pl]
    set angp [expr abs(asin($perp))]
    set pro [vecnorm [vecsub $rideal [vecscale $pl $perp]]]

    # compute decomposition angle
    set tang [expr cos($angp) * tan([expr 0.5 * $angl])]
    set angle [expr 2 * atan($tang)]

    # compute pivot point
    set hi [vecadd $bi [vecscale [veccross $pl $pro] [expr 0.5 * $cdis / $tang]]]

    # translate by effective rotation and reset
    set t [trans angle $com2 $hi $com1 $angle rad] 
    $all1 move $t
    set com3 [measure center $selid2_1 weight mass]
    set rmspro [measure rmsd $selid2_1 $selid2_2 weight mass]
    set t [measure fit $selid1_1 $selid1_2 weight mass]
    $all1 move $t
	
    # output 
    set pi [expr 2*asin(1.0)]
    set deg_angle [expr $angle * 180.0 / $pi] 
    set deg_angp [expr $angp * 180.0 / $pi] 
    set relerr [expr 100.0 * ($rmspro - $rmsid) / $cdis]
    puts "hinge> results:"
    puts "hinge> pivot point: $hi ."
    puts "hinge> effect. rot. axis: $pro (left-handed rot.)."
    puts "hinge> effect. rot. angle: $deg_angle degrees."
    puts "hinge> accuracy:"
    puts "hinge> rmsd (least squ.): $rmsid ."
    puts "hinge> rmsd (effect. rot.): $rmspro ."
    puts "hinge> relative error: $relerr %."
    puts "hinge> projection (deviation) angle: $deg_angp degrees."

    # draw hinge axis 
    set ax1 [vecadd $hi [vecscale $pro [expr 0.50 * $laxis]]]
    set ax2 [vecsub $hi [vecscale $pro [expr 0.42 * $laxis]]]
    set ax3 [vecsub $hi [vecscale $pro [expr 0.50 * $laxis]]]
    draw color [expr $domid2 % 17]
    draw cylinder $ax1 $ax2 radius $raxis
    draw cone $ax2 $ax3 radius [expr 2.5 * $raxis]
    draw sphere $ax1 radius $raxis
    draw cylinder $com1 $hi radius $raxis
    draw sphere $com1 radius $raxis
    draw cylinder $com2 $hi radius $raxis
    draw sphere $com2 radius $raxis
    return
}
