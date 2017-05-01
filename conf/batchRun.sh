#!/bin/bash

let index=0

for OPsim in 1 4 5 6 8 9 10 11
do

for OPrec in 1 4 5 6 8 9 10 11
do

    ./DarkSearch conf/config_C${OPsim}sim_C${OPrec}rec_200.dat > /dev/null

done
done
done
