# This measures the helix rotation between two trajectories
# according to the amide bond (C=O)

proc amide_rotation {hlx_beg hlx_end molid1 molid2 fid} {
	for {set k $hlx_beg} {$k <= $hlx_end} {incr k} {
		#find the C,O atoms
		set res_c1 [atomselect $molid1 "resid $k and type C"]
		set res_c2 [atomselect $molid2 "resid $k and type C"]
		set res_o1 [atomselect $molid1 "resid $k and type O"]
		set res_o2 [atomselect $molid2 "resid $k and type O"]
		
		#get their coordinates
		set c1 [lindex [$res_c1 get {x y z}] 0]
		set c2 [lindex [$res_c2 get {x y z}] 0]
		set o1 [lindex [$res_o1 get {x y z}] 0]
		set o2 [lindex [$res_o2 get {x y z}] 0]
		
		
		#create amide bond vectors
		set amide_vec1 [vecsub $o1 $c1]
		set amide_vec2 [vecsub $o2 $c2]
		
		#compute angle between the two amides
		set cos [vecdot [vecnorm $amide_vec1] [vecnorm $amide_vec2]]
		set rad_angle [expr {acos($cos)}]
		set deg_angle [expr  {$rad_angle*180/3.14159265359}]
		
		#print out
		puts $fid "[$res_c1 get resname] $deg_angle"
		
		##debug - check selections are OK
		#puts [$res_c1 get {resname type}]
		#puts [$res_c2 get {resname type}]
		#puts [$res_o1 get {resname type}]
		#puts [$res_o2 get {resname type}]

	}
}

proc rotation_by_helix_center {hlx_beg hlx_end molid1 molid2 fid} {
	for {set k $hlx_beg} {$k <= $hlx_end-4} {incr k} {
		set ca_sel1 [atomselect $molid1 "resid $k to [expr $k+3] and type CA"]
		set center_sel1 [measure center $ca_sel1]
		set third_ca1 [atomselect $molid1 "resid [expr $k+2] and type CA"]
		set p1 [lindex [$third_ca1 get {x y z}] 0]
		set vec1 [vecsub $p1 $center_sel1]
		
		set ca_sel2 [atomselect $molid2 "resid $k to [expr $k+3] and type CA"]
		set center_sel2 [measure center $ca_sel2]
		set third_ca2 [atomselect $molid2 "resid [expr $k+2] and type CA"]
		set p2 [lindex [$third_ca2 get {x y z}] 0]	
		set vec2 [vecsub $p2 $center_sel2]	
		
		set cos [vecdot [vecnorm $vec1] [vecnorm $vec2]]
		set rad_angle [expr {acos($cos)}]
		set deg_angle [expr  {$rad_angle*180/3.14159265359}]
		
		#print out
		puts $fid "[$third_ca1 get resname] $deg_angle"
		
	}
}
