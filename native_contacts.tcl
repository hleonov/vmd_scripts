source ~/vmd_scripts/utilities.tcl

proc compute_contacts {fr rad hx CAs seltxt ct_list} {
	upvar 1 $ct_list ct
	set counter 0
	$CAs frame $fr
	#puts "update CAs"
	foreach ca [$CAs list] {
       set sel [atomselect top "($seltxt) and (within $rad of (index $ca)) and not (same residue as (index $ca))" frame $fr]
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
					set ct($ca,$close_atom) $dist
				#	set ct($close_atom, $ca) $dist
					#puts "Contact: $ca ($rid1)\t$close_atom ($rid2)"
				}	
			}
	   }
	}
	#puts "passed compute_contacts"
	return $counter;
}

proc contact_list {fr rad CAs} {
	$CAs frame $fr
	set contact [measure contacts $rad $CAs $CAs]
	return [llength [lindex $contact 0]]
	
	
}

proc init_array {num arr} {
	upvar 1 $arr init
	for {set i 0} {$i < $num} {incr i} {
 	  for {set j 0} {$j < $num} {incr j} {
   		  set init($i,$j) 0;
  	  }
	}
	#puts "success in init!!"
}
# ============ Main ==================

check_for_help "vmd -dispdev text -e ~/vmd_scripts/native_contacts.tcl -args <gro> <out> <trj-type> <skip> <trj> <radius> <hx_dist> <native-pdb>"
#check_for_help "vmd -dispdev text -e ~/vmd_scripts/native_contacts.tcl -args <native>"

set pdb			[lindex $argv 0];
set out			[lindex $argv 1];
set trj_type	[lindex $argv 2];
set skip    	[lindex $argv 3]; 
set trj 		[lindex $argv 4];
set rad			[lindex $argv 5];
set hx_dist		[lindex $argv 6];
set native		[lindex $argv 7];
set seltxt 		[lreplace $argv 0 7];

#set seltxt "protein and noh and not backbone";
#set seltxt "protein and type CA";
#array unset nCAs_list;
set fid [open $out w];

#initialize array of native contacts 
set num 1000;
init_array $num nCAs_list;

#Load native (xray/nmr) structure
mol load pdb $native;

#compute native contacts
set nCAs [atomselect top "$seltxt"]
set num_cont [compute_contacts 0 $rad $hx_dist $nCAs $seltxt nCAs_list]

#array get nCAs_list
for {set i 0} {$i < $num} {incr i} {
   for {set j 0} {$j < $num} {incr j} {
     if {$nCAs_list($i,$j)>0} {
	     #puts "$i\t$j\t$nCAs_list($i,$j)";
		 
		 #eliminate duplicates
		 if {$nCAs_list($j,$i)>0} {
		 	set nCAs_list($j,$i) 0;
			set num_cont [incr num_cont -1]
		 }
	 }
   }
}
puts "native contacts number: $num_cont";
mol delete top 

#Load the reference structure and trajectory
mol load pdb $pdb;

if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}
set number_of_frames [molinfo top get numframes]
set CAs [atomselect top "$seltxt" frame 0];
##set close_list [atomselect top "protein and within $rad of ([$CAs text])" frame 0]

for {set fr 0} {$fr < $number_of_frames} {incr fr} {
	report $fr 50
	init_array $num temp_list;
	set temp_counter 0;
	set cont1 [compute_contacts $fr $rad $hx_dist $CAs $seltxt temp_list]
	for {set i 0} {$i < $num} {incr i} {
   		for {set j $i} {$j < $num} {incr j} {
     		if {($nCAs_list($i,$j)>0 && ($temp_list($i,$j)>0 || $temp_list($j,$i)>0)) ||
			    ($nCAs_list($j,$i)>0 && ($temp_list($i,$j)>0 || $temp_list($j,$i)>0))} {
	    		#puts $fid "$i\t$j\t$nCAs_list($i,$j)\t$temp_list($i,$j)";
				incr temp_counter;
	 		}
   		}
	}
	set perc [expr 100*$temp_counter/$num_cont];
	puts $fid "$fr\t$perc";
#	set cont2 [contact_list $fr $rad $CAs]
#	puts $fid "$fr\t $cont1 \t $cont2"
}

close $fid;
exit;
