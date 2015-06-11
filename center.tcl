# output the COM (center of mass) of the amantadine in each file (rc)

source /Users/hleonov/vmd_scripts/utilities.tcl
#check_for_help "vmdt -dispdev text -e ~/vmd_scripts/center.tcl -args [-h] <base> <out> <minZ> <maxZ> <dz> <use_end=1|0> <res> <k>"

set base [lindex $argv 0]
set out  [lindex $argv 1]
set min  [lindex $argv 2]
set max  [lindex $argv 3]
set dz   [lindex $argv 4]
set useend [lindex $argv 5]
set res		[lindex $argv 6]
set k 		[lindex $argv 7]

if {$useend == 1} {
	set end "_pull_k$k"
} else {
	set end ""
}
set fid [open $out w]
for {set i $min} {$i<=$max} {set i [expr $i+$dz]} {
	
	regsub {\.0} $i "" i2
	puts "i = $i\ti_sub=$i2"
	if {[file exists "$base\_$i2$end\.pdb"]} {
		mol load pdb "$base\_$i2$end\.pdb"
		set amn [atomselect top "resname $res"]
		puts "$res number: [$amn num]"
		set cent [measure center $amn]
		set cx [expr [lindex $cent 0]/10]
		set cy [expr [lindex $cent 1]/10]
		set cz [expr [lindex $cent 2]/10]
		puts $fid "$i2 $cx $cy $cz"
		mol delete top
	}
}
close $fid
exit;
