# Water connectivity:
# Writes for every frame if the pore has water connectivity in<-->out from given residue
# Usage example: vmd -dispdev text -e ~/vmd_scripts/water_connectivity.tcl -args system_from_PR2.gro system_from_MD.xtc outfile HISH skip radius

source /Users/hleonov/vmd_scripts/utilities.tcl
source /Users/hleonov/vmd_scripts/Tcl/bigdcd.tcl
source /Users/hleonov/vmd_scripts/Tcl/waitfor.tcl

check_for_help "Usage example: vmd -dispdev text -e ~/vmd_scripts/water_connectivity.tcl -args system_from_PR2.gro system_from_MD.trr connectivity HISH 10 3.5 trr"

set gro   [lindex $argv 0] ;	#starting structure 
set trj   [lindex $argv 1] ;	#trajectories 
set file  [lindex $argv 2] ;	#output 
set resname [lindex $argv 3]; 	#resname to search from
set skip_fr [lindex $argv 4];	#frame skipping when loading
set radius [lindex $argv 5]; 
set trjtype [lindex $argv 6];

# protein center +- 15A reaches bulk
set z_thresh 20	

mol load gro $gro
#animate read $trjtype $trj skip $skip_fr waitfor all
set fid [open $file w]


# -------------- Procs ------------------#
#get list of water OW within radius angstrom from sel
proc find_close_waters {sel fr} {
  global radius
  set near [atomselect top "type OW and within $radius of ([$sel text])" frame $fr]
#  set near [atomselect top "resname SOL HOH and within $radius of ([$sel text])" frame $fr]
  return [$near list]
}

#check if given atom z coordinate is in the bulk - plus side
proc beyond_plus_bulk {atom fr} {
	global z_thresh
	set tmp [atomselect top "index $atom" frame $fr]
	set z [$tmp get {z}]
	if {$z>=$z_thresh} {
		return 1
	}
	return 0
}
#check if given atom z coordinate is in the bulk - minus side
proc beyond_minus_bulk {atom fr} {
	global z_thresh
	set tmp [atomselect top "index $atom" frame $fr]
	set z [$tmp get {z}]
	if {$z<= [expr 0-$z_thresh]} {
		return 1
	}
	return 0
}

proc analyze_trj {fr} {
	global fid
	global resname
	center_by_selection "protein or resname NAC ACE" $fr
	#initial selection
	set sel [atomselect top "resname $resname and not backbone" frame $fr]
#	set sel [atomselect top "resname $resname" frame $fr]
	set global_water_list [] 
	report $fr 50
	set close_waters [find_close_waters $sel $fr]

	#set to 1 if reached bulk on each side
	set plus 0
	set minus 0
	while {([llength $close_waters] > 0) && ($plus<1 || $minus<1)} {
		set sel [atomselect top "index [$sel list]" frame $fr]
		set close_waters [find_close_waters $sel $fr]
		#remove duplicates
		foreach water $close_waters {
			set bulk 0
    		set i [lsearch $close_waters $water]
			set res [lsearch $global_water_list $water]
			#if current water is not already in global list - check bulk.
			if ($res<0) { 
				if {[beyond_plus_bulk $water $fr] == 1} {
					set plus 1
					set bulk 1
				} 
				if {[beyond_minus_bulk $water $fr] == 1} {
					set minus 1
					set bulk 1
				} 
			} 
			#if current water is already on global list or in  bulk - remove.
			if {($bulk == 1) || ($res>=0)} {
				set close_waters [lreplace $close_waters $i $i]
			}
    	}
		#add new water to global list
		foreach i $close_waters {
      		lappend global_water_list $i
    	}
		#update sel to be new water list
		#puts "$fr update to $global_water_list"
		set sel [atomselect top "index $global_water_list" frame $fr]
	}
	#echo "mol modselect 2 top \"index $global_water_list\""
	if {($plus==1) && ($minus==1) } {
		puts $fid "$fr\t1";
		#puts $fid "Frame $fr has a continuous water file YES"
		return 1;
	} else {
		puts $fid "$fr\t0";
		#puts $fid "Frame $fr - no continuous water file NO"
		return 0;
	}
}
proc end {} {
	exit;
}
# -------------------------- Main ------------------- #
waitfor $trj  end
bigdcd $trjtype analyze_trj $trj


###animate read $trjtype $trj skip $skip_fr waitfor all
#for {set fr 0} {$fr < [molinfo top get numframes]} {incr fr} {	
#	analyze_trj $fr
#}
#exit;
