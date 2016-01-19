#!/bin/bash

let index=0

for detector in "XE" "AR" "GE" "SI" "FL"
do

for OP in C1 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15
do

let index+=1
    if [ ! -f ../../Projects/nuFloor/bestFit_NRops/data/DS_sim_${detector}_${OP}_n.dat ]; then
        
        echo "running $index"
        ./DarkSearch conf/config${index}.dat

        mv results/DS_sim.dat ../../Projects/nuFloor/bestFit_NRops/data/DS_sim_${detector}_${OP}_n.dat
        mv results/DS_stats.dat ../../Projects/nuFloor/bestFit_NRops/data/DS_stats_${detector}_${OP}_n.dat
    fi

done
done
