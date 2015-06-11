source ~/vmd_scripts/utilities.tcl

proc compute_contacts {fr rad hx CAs close_list} {
	set counter 0
	$CAs frame $fr
	#puts "update CAs"
	foreach ca [$CAs list] {
       set sel [atomselect top "(protein and type CA) and (within $rad of (index $ca)) and not (same residue as (index $ca))" frame $fr]
	   set ca_sel [atomselect top "index $ca" frame $fr]
	 #  puts "frame $fr: CA is $ca"
	   foreach close_atom [$sel list] {
	   		set close_sel [atomselect top "index $close_atom" frame $fr]
		#	puts [$sel num]
		#	puts "frame $fr: close atom is $close_atom"
			set rid1 [$ca_sel get resid]
			set rid2 [$close_sel get resid]
			if {[expr abs($rid1-$rid2)] >= $hx} {
				set coord1 [lindex [$ca_sel get {x y z}] 0];
				set coord2 [lindex [$close_sel get {x y z}] 0];
				set dist [veclength [vecsub $coord1 $coord2]];
	
			#	puts "rid1= $rid1 \t rid2 = $rid2 \t Distance = $dist"
				if {$dist<$rad} {
					incr counter
				#	puts "counter: $counter"
				}	
			}
	   }
	}
	return [expr $counter/2]
}

proc contact_list {fr rad CAs} {
	$CAs frame $fr
	set contact [measure contacts $rad $CAs $CAs]
	return [llength [lindex $contact 0]]
	
	
}
# ============ Main ==================

check_for_help "vmd -dispdev text -e ~/vmd_scripts/contacts.tcl -args <gro> <out> <trj-type> <skip> <trj> <radius> <hx_dist>"

set pdb			[lindex $argv 0];
set out			[lindex $argv 1];
set trj_type	[lindex $argv 2];
set skip    	[lindex $argv 3]; 
set trj 		[lindex $argv 4];
set rad			[lindex $argv 5];
set hx_dist		[lindex $argv 6];

set CAs 0
set close_list 0

mol load pdb $pdb;
set fid [open $out w];
if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}
set number_of_frames [molinfo top get numframes]
set CAs [atomselect top "protein and type CA" frame 0];
#set close_list [atomselect top "protein and within $rad of ([$CAs text])" frame 0]

for {set fr 0} {$fr < $number_of_frames} {incr fr} {
	report $fr 50
	set cont1 [compute_contacts $fr $rad $hx_dist $CAs $close_list]
	set cont2 [contact_list $fr $rad $CAs]
	puts $fid "$fr\t $cont1 \t $cont2"
}

close $fid;
exit;
