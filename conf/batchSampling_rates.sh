#!/bin/bash

OPsVal=(1.78e-4 0.192 6.05 52.1 0.353 5.27 2.06 .00606)
detRatio=(.195 .183 .270 .786 .186 .187 .446 .287)
index=0

for OPsim in 1 4 5 6 8 9 10 11
do

for OPrec in 1 4 5 6 8 9 10 11
do

cat > config_C${OPsim}sim_C${OPrec}rec_20a.dat << EOF
//Dark Search running options
0                  // running mode: 0 simulation, 1 Print out dR/dE
./results/DS_C${OPsim}sim_C${OPrec}rec_20a_      // root of output file names, can include directories, make unique to run 2 instances of DS simultaneously
//Multinest Sampling Parameters (0 off, 1 on)
1                  // do mode separation?
0                  // run in constant efficiency mode?
4000                // number of live points
1                  // need feedback to standard output (terminal)?
0                  // resume from a previous job? (only available if asimov data being used)
0.5                // set the required efficiency. 0.8 and 0.3 are recommended for parameter estimation & evidence evalutaion respectively.
0.1                // tol, defines the stopping criteria, 0.5 gives good enough accuracy
-1.e90             // all the modes with logZ < Ztol are ignored
0                  // Use a binless likelihood (0=no/1=yes)
//Include coherent neutrino scattering background
0
//Allow isospin violation in reconstruction? (0=no/1=yes)
0
//Search range of parameters
//par  | min/mean   | max/error   | prior  (log/linear/gaussian/none(mean used))
Mx       1            1e3             log           // WIMP mass (GeV)  
spin     0.5                          none
C1       1e-5         1e2             none                  
C2       1e-5         1e2             none
C3       1e-5         1e2             none
C4       1e-5         1e2             none
C5       1e-5         1e2             none
C6       1e-5         1e2             none
C7       1e-5         1e2             none
C8       1e-5         1e2             none
C9       1e-5         1e2             none
C10       1e-5         1e2             none
C11       1e-5         1e2             none
C12       1e-5         1e2             none
C13       1e-5         1e2             none
C14       1e-5         1e2             none
C15       1e-5         1e2             none
rho      0.3          0.1           gaussian       // Local dark matter density (Gev/cm^3)
v0       220          20            gaussian       // Galactic rotation velocity (km/s)
vesc     544          40            gaussian       // Galactic escape velocity from earth (km/s)
vSp      10           0             none       // Magnitude of Solar peculiar velocity w.r.t local standard of rest
vEp      0            0             none       // Magnitude of Earth peculiar velocity w.r.t the Sun
// Detector  |  Exposure (Tonne*years)
#XENON         .1
#GERMANIUM     ${detRatio[$index]}
//Simulated WIMP parameters
Mx       100
spin     0.5        
//operator  protons   neutrons
C1          0          0
C2          0          0
C3          0          0
C4          0          0
C5          0          0
C6          0          0
C7          0          0
C8          0          0
C9          0          0
C10          0          0
C11          0          0
C12          0          0
C13          0          0
C14          0          0
C15          0          0
rho     0.3          
v0      220          
vesc    544          
vSp      10           
vEp       0 
//Simulate Asimov data set/random Monte-Carlo (0/1)
0
EOF

val=${OPsVal[$index]}
sed -i "s|C$OPrec       1e-5         1e2             none|C$OPrec       1e-5         1e2            log|g" config_C${OPsim}sim_C${OPrec}rec_20a.dat 
sed -i "s|C$OPsim          0          0|C$OPsim    ${OPsVal[$index]}   ${OPsVal[$index]}|g" config_C${OPsim}sim_C${OPrec}rec_20a.dat 

done
let index+=1
done
