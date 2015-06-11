source /Users/hleonov/vmd_scripts/utilities.tcl
check_for_help "vmdt -dispdev text -e ./align_drg_z.tcl -args <drug.pdb> <atom1> <atom2> <1|0 (down|up)>"
set drg_file [lindex $argv 0]
set a1 		 [lindex $argv 1]
set a2		 [lindex $argv 2]
set down	 [lindex $argv 3]; 		#if down=1:  set A [orient $p_drg [lindex $I 2] {0 0 -1}]
									#else   	 set A [orient $p_drg [lindex $I 2] {0 0 1}]
#set drg_file "Riman_prodrgB.pdb"
#set a1 CAJ
#set a2 CAK

mol load pdb $drg_file
set drug [atomselect top "all"]
$drug num

#move system to be centered at 0,0,0
$drug moveby [vecinv [measure center $drug]]

set p_drg [atomselect top "type $a1 $a2"]
$p_drg num

#align drug vector with z-axis
package require Orient
namespace import Orient::orient

set I [draw principalaxes $p_drg]
if {$down==1} {
	set A [orient $p_drg [lindex $I 2] {0 0 1}]	
} else {
   set A [orient $p_drg [lindex $I 2] {0 0 -1}]
}
$drug move $A
set I [draw principalaxes $p_drg]

namespace forget Orient::orient
if {$down==1} {
	$drug writepdb "rot_down_$drg_file"
} else {
	$drug writepdb "rot_up_$drg_file"
}
exit
#mol delete top

