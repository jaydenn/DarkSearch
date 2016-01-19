#!/bin/bash

let index=0

for detector in "XENON_high" "GERMANIUM_high"
do

for OP in C1 C6 C10
do

for N in n p
do

    ./DarkSearch conf/config_${detector:0:2}_${OP}_${N}.dat > /dev/null

    mv results/DS_${detector}_dRdE.dat ../../Projects/nuFloor/bestFit_NRops/results_atm/DS_rates_${detector:0:2}_${OP}_${N}.dat

done
done
done
