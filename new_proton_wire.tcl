# Water connectivity:
# Writes for every frame if the pore has water connectivity in<-->out from given residue
# Usage example: vmd -dispdev text -e ~/vmd_scripts/water_connectivity.tcl -args system_from_PR2.gro system_from_MD.xtc outfile HISH skip radius

source /Users/hleonov/vmd_scripts/utilities.tcl

check_for_help "Usage example: vmd -dispdev text -e ~/vmd_scripts/new_proton_wire.tcl -args <in.pdb> <output> <trj.xtc> 10 xtc 3.5 1 resname HIS"

set pdb   [lindex $argv 0] ;	#starting structure 
set file  [lindex $argv 1] ;	#output 
set trj   [lindex $argv 2] ;	#trajectories 
set skip  [lindex $argv 3];		#frame skipping when loading
set trjtype [lindex $argv 4];	#xtc or trr
set cutoff  [lindex $argv 5]; 	#cutoff radius
set use_n	[lindex $argv 6];  ; #use N atoms
set res_sel [lreplace $argv 0 6]; 	#residue expression to search from

set bulk_thresh 5	; # reached bulk when more than $bulk_thresh water molecules are around 
mol load pdb $pdb
if {($skip>1)} {
	animate read $trjtype $trj skip $skip waitfor all
} else {
	animate read $trjtype $trj waitfor all
}
set fid [open $file w]


# -------------- Procs ------------------#

proc check_bulk {fr front_list start_Z} {
	global bulk_thresh cutoff
	foreach atom $front_list {
		set neig_sel [atomselect top "type OW and (within $cutoff of index $atom) and not index $atom" frame $fr]
		if {[$neig_sel num] >= $bulk_thresh} {
			#up or down?
			set tmp [atomselect top "index $atom" frame $fr]
			set z [$tmp get {z}]
			if {$z < $start_Z} {
				set minus 1
			} else { 
				set plus 1 
			}
			
			#remove
			set elm_i [lsearch $front_list $atom]
			#puts "remove $atom, index $elm_i"
			if {$elm_i > -1} {
				set front_list [lreplace $front_list $elm_i $elm_i]
			}
		}
	}
	return $front_list
}

proc get_h_bond_list {frame} {
  global res_sel cutoff use_n argv plus minus
  set start_atoms [atomselect top "protein and $res_sel and sidechain and type \"O.*\" \"N.*\"" frame $frame]
  set start_Z [lindex [measure center $start_atoms] 2]
  set current [$start_atoms list]; #get indices list of starting  protein atoms
  set cummulative $current
#	set i 0
#	set repeats 1
  while {($plus == 0) && ($minus == 0)} {
  
  
  #while {$i < $repeats} {}
    if {$use_n} {
      set front [atomselect top "(type \"O.*\" \"N.*\" and within $cutoff of (index $current)) and not (index $cummulative $current)"  frame $frame]    
    } else {
      set front [atomselect top "(type \"O.*\" and within $cutoff of (index $current)) and not (index $cummulative $current)"  frame $frame]    
    }
    set cummulative [concat $cummulative $current];
    set current [$front list];
#	puts "before: [llength $current]"

    if {$current == 0 || $current == "" } {
      break;
    }
#	incr i

	set current [check_bulk $frame $current $start_Z]
#	puts "after:  [llength $current]"
  }
  return $cummulative;
}

##--- main -- ###
for {set i 1} {$i < [molinfo top get numframes]} {incr i} {
	if { $frame % 50 == 0 } {
	    puts "analyzing frame $i"
	}
    puts $fid [get_h_bond_list $i];
}
