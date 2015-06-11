set idir /data/desrad-nb-o/jensenm/Aqp0/2B6Orun/desmondNPT_p40614/NPT_Berendsen_early_equil_2fs/
set xdir /data/desrad-nb-o/jensenm/Aqp0/2B6Orun/desmondNPT_p40614/NPT_Berendsen_early_equil_2fs/xtc
set rdir /data/desrad-nb-o/jensenm/Aqp0/2B6Orun/desmondNPT_p40614/NPT_Berendsen_early_equil_2fs/analysis/HOLE
set wdir ${rdir}
set bindir /data/desrad-nb-o/jensenm/bin
set f [expr 1.0 / 3.0]

mol load pdb ${idir}/pdb/Aqp0_POPE_frame_1_cleaned_w_segnames.pdb
set all  [atomselect top "protein"]
set ref  [atomselect top "protein and name CA" frame 0]
set cmp  [atomselect top "protein and name CA" ]

set off 833 ; # 5 ns 
set skip 5000 ; # in ps
set xtclist [list "0-157272.xtc"]

set xtc [lindex $xtclist 0]

set ncyc [expr int(floor( (157272-$skip) / $skip))]


for {set k 0} {$k < [expr $ncyc -1]} {incr k} {

    set beg [expr $off + $k*$off]
    set fin [expr $beg + $off]

    animate read xtc ${xdir}/$xtc beg $beg skip 5 end $fin waitfor all

    set pro1 [atomselect top "segname PRO1 and not hydrogen"]
    set pro2 [atomselect top "segname PRO2 and not hydrogen"]
    set pro3 [atomselect top "segname PRO3 and not hydrogen"]
    set pro4 [atomselect top "segname PRO4 and not hydrogen"]

    for {set i 0} {$i < [molinfo top get numframes]} {incr i} {

	$all frame $i
	$all move [measure fit $cmp $ref]

        $pro1 frame $i
        $pro1 writepdb p1_${i}.pdb
        $pro2 frame $i
        $pro2 writepdb p2_${i}.pdb
        $pro3 frame $i
        $pro3 writepdb p3_${i}.pdb
        $pro4 frame $i
        $pro4 writepdb p4_${i}.pdb

        for {set m 1 } {$m < 5} {incr m} {

            set fid [open hole.inp w+]
            puts "DCD: $k; grepping data from p${m}_${i}.pdb"

            set r1 [exec grep PHE ${rdir}/p${m}_${i}.pdb | grep CG  | grep " 48 "  | cut -c30-54]
            set r2 [exec grep HIS ${rdir}/p${m}_${i}.pdb | grep NE2 | grep " 172 " | cut -c30-54]
            set r3 [exec grep ARG ${rdir}/p${m}_${i}.pdb | grep CZ  | grep " 187 " | cut -c30-54]

            set cpoint {}

            for {set n 0} {$n < 3} {incr n} {
                set c${n} "[expr $f * ([lindex $r1 $n] + [lindex $r2 $n] + [lindex $r3 $n] )]"
            }

            lappend cpoint $c0 $c1 $c2

            puts $cpoint
            puts $fid "coord  p${m}_${i}.pdb"
            puts $fid "radius simple.rad"
            puts $fid "cvect  0 0 1"
            puts $fid "sample 0.1"
            puts $fid "endrad 5"
            puts $fid "sphpdb ${k}_p${m}_${i}.hole.sph"
            puts $fid "cpoint $cpoint"
            flush $fid
            close $fid
            exec $bindir/hole < hole.inp 
            exec rm -f hole.inp
            exec rm -f p${m}_${i}.pdb
        }
    }
    animate delete all
}

