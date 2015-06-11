source ~/vmd_scripts/utilities.tcl
source ~/vmd_scripts/rmsdalign.tcl

# Reads a gro file and following trajectories, and computes the Z distribution of each residue along the pore.
# The average is given on all 4 helices together. 


check_for_help "vmd -dispdev text -e ~/vmd_scripts/z_aa.tcl -args <pdb> <trj-type> <skip> <trj> <chain_length> <0|1|2=no|incl|only backbone>"

set pdb			[lindex $argv 0];
set trj_type	[lindex $argv 1];
set skip    	[lindex $argv 2]; 
set trj 		[lindex $argv 3];
set chain_length	[lindex $argv 4];
set bb			[lindex $argv 5];

mol load pdb $pdb;

# set the selection arrays and file handle arrays

for {set i 1} {$i <= $chain_length} {incr i} {
  if {$bb == 0} { ; # only sidechains
	  set selections($i) [atomselect top "resid $i [expr $i + (1 * $chain_length)] [expr $i + (2 * $chain_length)] [expr $i + (3 * $chain_length)] and not backbone"]
  } elseif {$bb == 1} {	; # sidechain+backbone
      set selections($i) [atomselect top "resid $i [expr $i + (1 * $chain_length)] [expr $i + (2 * $chain_length)] [expr $i + (3 * $chain_length)]"]
  } elseif {$bb == 2} {	; # only backbone
	set selections($i) [atomselect top "resid $i [expr $i + (1 * $chain_length)] [expr $i + (2 * $chain_length)] [expr $i + (3 * $chain_length)] and backbone"]
  }
  set name [lindex [$selections($i) get resname] 0]
  set out($i) [open "res.${name}_${i}.dat" w]
}

# now read the trajectory file

if {($skip>1)} {
	animate read $trj_type $trj skip $skip waitfor all
} else {
	animate read $trj_type $trj waitfor all
}

#align according to CA
#set CA_ref [atomselect top "protein and type CA" frame 0]
#set CA_sel [atomselect top "protein and type CA" frame 0]
#set all_system [atomselect top "all"]
#align $CA_sel $CA_ref $all_system top

#get Z of aminoacids
for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
    for {set j 1} {$j <= $chain_length} {incr j} {
      $selections($j) frame $i
      set z [lindex [measure center $selections($j)] 2]
      puts $out($j) "$z"
	  flush $out($j); # not really needed
    }
  }

for {set i 1} {$i <= $chain_length} {incr i} {
  close $out($i)
}

unset selections
exit


