#!/usr/local/bin/bash

mono=1
nm=4
nruns=12

while [ $mono -le $nm ]
do
	run=0
	while [ $run -le $nruns ]
	do
 		hol=0_p${mono}_${run}.hole.sph
                if [ -e ${hol} ]; then
                	echo "${hol} present, extracting radii"
			rm -f ${mono}_${run}.dat
			cut -c48-60 $hol > tmp/${mono}_${run}.dat
		else	
			echo "HOLE output not present"
		fi
		run=$((run+1))
	done
	mono=$((mono+1))
done
