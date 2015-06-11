# run like: vmdt -e check_close.tcl -args ABZ_solvated.pdb 3
 set pdb  [lindex $argv 0]
set threshold   [lindex $argv 1]
mol load pdb $pdb
set sel [atomselect top "(resname ARG and (name NE NH1 NH2)) or (resname LYS and name NZ)"];
foreach i [$sel list] {
  set tmpi [atomselect top "index $i"]
  set veci [lindex [$tmpi get {x y z}] 0]
  set resi [$tmpi get resid]
  set resiname [$tmpi get resname]
  set contacts [lindex [measure contacts $threshold $tmpi $sel] 1]
  if {[llength $contacts] > 0} {
    foreach j $contacts {
      set tmpj [atomselect top "index $j"]
      set vecj [lindex [$tmpj get {x y z}] 0]
      set resj [$tmpj get resid]
      set resjname [$tmpj get resname]
      if {$resi != $resj} {
        set d [vecdist $veci $vecj]
        puts "Residue: ${resiname}-${resi} is within $d of Residue: ${resjname}-${resj}"
      }
    }
  }
}
set sel [atomselect top "(resname GLU and (name OE1 OE2)) or (resname ASP and (name OD1 OD2))"];
foreach i [$sel list] {
  set tmpi [atomselect top "index $i"]
  set veci [lindex [$tmpi get {x y z}] 0]
  set resi [$tmpi get resid]
  set resiname [$tmpi get resname]
  set contacts [lindex [measure contacts $threshold $tmpi $sel] 1]
  if {[llength $contacts] > 0} {
    foreach j $contacts {
      set tmpj [atomselect top "index $j"]
      set vecj [lindex [$tmpj get {x y z}] 0]
      set resj [$tmpj get resid]
      set resjname [$tmpj get resname]
      if {$resi != $resj} {
        set d [vecdist $veci $vecj]
        puts "Residue: ${resiname}-${resi} is within $d of Residue: ${resjname}-${resj}"
      }
    }
  }
}
exit
