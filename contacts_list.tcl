# Create the elastic network model list
# as an input to the system.gmx file.

# Initialize psf and coordinates
mol load pdb basic.pdb

# Select the type of atoms on which the network will be built
#set ENMatom "name N C O CA" 
set ENMatom "protein and resid 150 to 240 and noh"
# Alternatives include noh, or main chain atoms...

# The network cutoff, network elastic constant
set cutoff 6.0
set ENMk   1000.0

# functype index and output file
# Set this variable by the output of grep ntypes system.gmx
# set indexFunctype 433
set indexFunctype [exec grep ntypes system.gmx | cut -f2 -d=]
set fidFunctype [open "functype.add" w+]

# bond index and output file
# Set this variable to the current number of bonds --- the output of:
#set indexBond 11841
set indexBond [expr [exec grep -A1 " Bond:" system.gmx | grep nr | cut -f2 -d: ] / 3]
set fidBond [open "bond.add" w+]

# Select the atoms on which the network will be built
set sel [atomselect top $ENMatom]
set listSel [$sel list]
set nsel [$sel num]

# Get the global coordinate matrix
set all [atomselect top "all"]
set cor [$all get {x y z}]
$all delete

# Get all the contacts
# first delete all the bonds
set elist {{}}
for {set i 0} {$i < $nsel} {incr i} {lappend elist {}}
$sel setbonds $elist
# use vmd's build in measure contacts facility:
set contact [measure contacts $cutoff $sel $sel]

# List all the pairs:
foreach i [lindex $contact 0] j [lindex $contact 1] {
    set r1 [lindex $cor $i]
puts "$r1"
    set r2 [lindex $cor $j]
puts "$r2"
    set dist [expr [vecdist $r1 $r2] / 10.0]
    puts $dist
    puts $fidBond     "               $indexBond type=$indexFunctype (BONDS) $i $j"
    puts $fidFunctype "            functype\[$indexFunctype\]=BONDS, b0A= $dist, cbA= $ENMk, b0B= $dist, cbB= $ENMk"
    incr indexFunctype
    incr indexBond
    flush $fidBond
    flush $fidFunctype 
}
$sel delete
close $fidFunctype
close $fidBond

# Now a major hack to create a new system.gmx file within tcl.
# Please manually verify the result as certain assumptions about
# the sequence of tags in system.gmx may not hold true in your system.

set fidFunctype [open "functype.add" r]
set fidBond [open "bond.add" r]

set sysgmx [open "system.gmx" r]
set sysenm [open "system-enm.gmx" w+]
set afterNtypes 0
set done 0
while {1} {
    if {$done==1} {
	break
    }

    gets $sysgmx line
    if {[eof $sysgmx]} {
	close $sysgmx
	break
    }
    
    if { [regexp {ntypes} $line] } {
	puts $sysenm "         ntypes=$indexFunctype"
	set afterNtypes 1
    } elseif { $afterNtypes==1 && [regexp {G96Bond:} $line] } {
	puts $sysenm [read -nonewline $fidBond]
	close $fidBond
	puts $sysenm $line
	# Past the remaining part of the file and clean up.
	while {1} {
	    gets $sysgmx line
	    puts $sysenm $line
	    if {[eof $sysgmx]} {
		close $sysgmx
		set done 1
		break
	    }
	}
    } elseif { $afterNtypes==1 && [regexp {Bond:} $line] } {
	puts $sysenm [read -nonewline $fidFunctype]
	close $fidFunctype
	puts $sysenm $line
	gets $sysgmx line
	puts $sysenm "            nr: [expr 3*($indexBond)]"
	gets $sysgmx line
	puts $sysenm "            multinr\[division over processors\]: [expr 3*($indexBond)]"
    } else {
	puts $sysenm $line
    }
}

close $sysenm

exit







