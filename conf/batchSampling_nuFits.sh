#!/bin/bash

let index=0

for detector in "#XENON_low  .22" "#ARGON_low  .3" "#GERMANIUM_low  0.15" "#SILICON_low  0.55" "#FLOURINE_low 1.8"
do

for OP in C1 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15
do

let index+=1
cat > config$index.dat << EOF
//Darth Sidious running options
0                  // running mode: 0 simulation, 1 Print out dR/dE
./results/DS_      // root of output file names, can include directories, make unique to run 2 instances of DS simultaneously
//Multinest Sampling Parameters (0 off, 1 on)
1                  // do mode separation?
0                  // run in constant efficiency mode?
1000                // number of live points
1                  // need feedback to standard output (terminal)?
0                  // resume from a previous job? (only available if asimov data being used)
0.8                // set the required efficiency. 0.8 and 0.3 are recommended for parameter estimation & evidence evalutaion respectively.
0.5                // tol, defines the stopping criteria, 0.5 gives good enough accuracy
-1.e90             // all the modes with logZ < Ztol are ignored
0                  // Use a binless likelihood (0=no/1=yes)
//Include coherent neutrino scattering background
1
//Search range of parameters
//par  | min/mean   | max/error   | prior  (log/linear/gaussian/none(mean used))
Mx       1            20              linear           // WIMP mass (GeV)  
C1       1e-8         1e3             none                  
C2       1e-8         1e3             none
C3       1e-8         1e3             none
C4       1e-8         1e3             none
C5       1e-8         1e3             none
C6       1e-8         1e3             none
C7       1e-8         1e3             none
C8       1e-8         1e3             none
C9       1e-8         1e3             none
C10       1e-8         1e3             none
C11       1e-8         1e3             none
C12       1e-8         1e3             none
C13       1e-8         1e3             none
C14       1e-8         1e3             none
C15       1e-8         1e3             none
rho      0.3          0.1           none //gaussian       // Local dark matter density (Gev/cm^3)
v0       220          20            none //gaussian       // Galactic rotation velocity (km/s)
vesc     544          40            none //gaussian       // Galactic escape velocity from earth (km/s)
vSp      10           0             none       // Magnitude of Solar peculiar velocity w.r.t local standard of rest
vEp      0            0             none       // Magnitude of Earth peculiar velocity w.r.t the Sun
// Detector  |  Exposure (Tonne*years)
${detector}     
//Simulated WIMP parameters
Mx      6        // WIMP mass [GeV]
spin    0.5        // WIMP spin
C1      4.5e-40                           
C2       0
C3       0           
C4       0           
C5       0           
C6       0           
C7       0           
C8       0           
C9       0           
C10      0           
C11      0           
C12      0           
C13      0           
C14      0           
C15      0           
//Simulate Asimov data set/random Monte-Carlo (0/1)
0
EOF

#eval "sed -i -e 's/"$OP"       1e-8         1e3             none/"$OP"       1e-8         1e3             log/g' $1"
sed -i "s|$OP       1e-8         1e3             none|$OP       1e-8         1e3             log|g" config$index.dat 

done
done
