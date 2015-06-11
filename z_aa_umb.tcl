# read in the molecule
set k [lindex $argv 0]
mol load gro center_0.0_pull_k$k.gro

# set the selection arrays and file handle arrays
set chain_length 27
for {set i 1} {$i <= $chain_length} {incr i} {
  set selections($i) [atomselect top "resid $i [expr $i + (1 * $chain_length)] [expr $i + (2 * $chain_length)] [expr $i + (3 * $chain_length)] and not backbone"]
  set name [lindex [$selections($i) get resname] 0]
  set out($i) [open "res.${name}_${i}.dat" w]
}

# now read in all the XTC files one by one

for {set x -14.5} {$x <= 14.5} {set x [expr $x + 0.5]} {
  set xtc "center_${x}_pull_k$k.xtc"
  puts $xtc
  animate read xtc $xtc waitfor all
  for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
    for {set j 1} {$j <= $chain_length} {incr j} {
      $selections($j) frame $i
      set z [lindex [measure center $selections($j)] 2]
      puts $out($j) "$z"
	  flush $out($j); # not really needed
    }
  }
  animate delete beg 1
}
for {set i 1} {$i <= $chain_length} {incr i} {
  close $out($i)
}

unset selections
exit
