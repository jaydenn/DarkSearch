#!/bin/bash

let index=0

for detector in "XENON_low" "GERMANIUM_low" "SILICON_low" "FLOURINE_low"
do

for OP in C1 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15
do

for N in n p
do

    ./DarkSearch conf/config_${detector:0:2}_${OP}_${N}.dat > /dev/null

    mv results/DS_${detector}_dRdE.dat ../../Projects/nuFloor/bestFit_NRops/results/DS_rates_${detector:0:2}_${OP}_${N}.dat

done
done
done
