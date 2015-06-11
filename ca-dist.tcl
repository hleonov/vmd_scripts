# given the selection, get the CAs and draw a distance plot


#VMD  --- start of VMD description block
#Name:
#  CA distance plot
#Synopsis:
#  Given a selection, plot the CA-CA distance plot
#Version:
#  1.1
#Uses VMD version:
#  1.5 
#Ease of use:
#  2
#Procedures:
# <li> ca_dist selection -- makes a new graphics molecule which is
#   the CA-CA distance plot of the given selection
#Description:
#   Given a selection, this script finds the CA atoms and constructs a
# matrix whose elements are the distances from the CA of residue(i) to
# that on residue(j).  A new graphics molecule is created with the name
# "CA_distance(molid)", where molid is the index of the molecule used
# by the selection.  The distance matrix is plotted by residue number
# in the order they appeared in the PDB (or PSF, ...) file, but only
# the selected residues are shown.  The colors are determined by the
# color scale and may be changed accordingly.
#Example output:
# <a href="br.dist.gif"><img src="br.dist.small.gif">image</a>
# of bacteriorhodopsin with the CA-CA  on the bottom
#Files: 
# <a href="ca-dist.vmd">ca-dist.vmd</a>
#Authors: 
# Andrew Dalke &lt;dalke@ks.uiuc.edu&gt;
# Justin Gullingsrud &lt;justin@ks.uiuc.edu&gt;
#\VMD  --- end of block

# Input: a selection
# Does: finds the CAs in the selection then computes and draws the
#   CA-CA distance grid with colors based on the color scale
# Returns: the id of the new graphics molecule containing the grid

proc mymax { a b } {
  if {[expr $a > $b]} {
    return $a
  } else {
    return $b
  }
}

proc ca_distance {main_sel} {
    # get the CA atoms from the selection
    set mol [$main_sel molindex]
    set seltext [$main_sel text]
    set sel [atomselect $mol "($seltext) and name CA"]

    # find distances between each pair
    set coords [$sel get {x y z}]
    set max 0
    set list2 [$sel list]

    foreach atom1 $coords id1 [$sel list] {
	foreach atom2 $coords id2 $list2 {
	    set dist($id1,$id2) [veclength [vecsub $atom2 $atom1]]
	    set dist($id2,$id1) $dist($id1,$id2)
	    set max [mymax $max $dist($id1,$id2)]
	}
	lvarpop list2
	lvarpop coords
    }


    # draw the pretty graphic
    puts "Distances calculated, now drawing the distance map ..."
    mol load graphics "CA_distance($mol)"
    set gmol [molinfo top]
    # turn material characteristics off
    graphics $gmol materials off
    # i1 and j1 are i+1 and j+1; this speeds up construction of x{01}{01}
    set i 0
    set i1 1
    # preset the scaling factor based on the number of available map colors
    set mapcolors [expr [colorinfo max] - [colorinfo num]]
    set numcolors [colorinfo num]
    set scale [expr ($mapcolors - 0.05) / ($max + 0.)]
    set list2 [$sel list]
    foreach id1 [$sel list] {
	set j 0
	set j1 1
	set x00 [list $i $j 0]
	set x10 [list $i1 $j 0]
	set x11 [list $i1 $j1 0]
	set x01 [list $i $j1 0]
	foreach id2 $list2 {
	    set col [expr int($scale * $dist($id1,$id2)) + $numcolors]
	    graphics $gmol color $col
	    graphics $gmol triangle $x00 $x10 $x11
	    graphics $gmol triangle $x00 $x11 $x01
	    incr j
	    incr j1
	    set x00 $x01
	    set x10 $x11
	    set x11 [list $i1 $j1 0]
	    set x01 [list $i $j1 0]
	}
	incr i
	incr i1
    }
    # put some numbers down to give an idea of what is where
    set resids [$sel get resid]
    set num [llength $resids]
    set start [lindex $resids 0]
    set end [lindex $resids [expr $num - 1]]
    graphics $gmol color yellow
    graphics $gmol text "0 0 0" "$start,$start"
    graphics $gmol text "$num 0 0" "$end,$start"
    graphics $gmol text "0 $num 0" "$start,$end"
    graphics $gmol text "$num $num 0" "$end,$end"
    return $gmol
}
