# Move bilayer to be centered at {0,0,0}
# Align pore with Z-axis
# Move protein to the center of the bilayer
# was saved as  system_pbc_tr_nohole.pdb (for no whole dmpc) system_tr_nohole.pdb
# then in visual vmd - select the non-cliding atoms using the selection txt in hole_sel.txt
# and save as system_pbc_hole.pdb system_hole.pdb

source ~/vmd_scripts/utilities.tcl
check_for_help "vmdt -e ~/vmd_scripts/align.tcl -args <bilayer-pdb> <prot-pdb>"
set bilayer_pdb [lindex $argv 0]
set prot_pdb [lindex $argv 1]

mol load pdb $bilayer_pdb
set bilayer [atomselect top "all"]

#protein (+HOH?)
mol load pdb $prot_pdb
set protein [atomselect top "all"]
set all_prot [atomselect top "all"]

$bilayer num
$protein num

#move system to be centered at 0,0,0
$bilayer moveby [vecinv [measure center $bilayer]]

#align pore with z-axis
package require Orient
namespace import Orient::orient

set I [draw principalaxes $protein]
set A [orient $protein [lindex $I 2] {0 0 1}]
$all_prot move $A
set I [draw principalaxes $protein]
set A [orient $protein [lindex $I 1] {0 1 0}]
$all_prot move $A
set I [draw principalaxes $protein]

namespace forget Orient::orient
$all_prot moveby [vecinv [measure center $protein]]

set origin [measure center $bilayer]
set prot_vec [measure center $protein]
set dist_vec [vecsub $origin $prot_vec]

$all_prot moveby $dist_vec

$all_prot writepdb "prot_align.pdb"
$bilayer writepdb "bilayer_align.pdb"

exit
#mol delete top
#mol delete top
