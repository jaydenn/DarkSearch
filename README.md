# Dark Search

June 2015

This code contains functions for the calculation of direct detection rates with full nuclear form factors (calculated by DMformfactor). It is designed to simulate a WIMP signal and reconstruct WIMP model parameters from the simulated signal.


## Installation
requires:
- Multinest         
 -- nested sampling algorithm for posterior sampling (requires lapack) 
 -- https://ccpforge.cse.rl.ac.uk/gf/project/multinest/
- GNU Science Lib   
 -- probability functions, integration and random number generation.
 -- http://www.gnu.org/software/gsl/

On debian based systems type:
```
sudo apt-get install libgsl0-dev liblapack-dev gfortran
```

optional:
openmpi - needed for parallelization of sampling 
```
sudo apt-get install openmpi-bin libopenmpi-dev 
```
To build just type: 
```
make
```

To build with mpi: 
1. modify the makefile to use mpif90 and uncomment the extra mpi flags
2. ensure Multinest was compiled with mpi support
2. type 'make' (run 'make clean' first if previously built without mpi) 

After installation it would be wise to type:
```
git update-index --assume-unchanged Makefile config.dat
```

this will prevent these files from being overwritten on future pulls from the server.

## config.dat
All of the running/sampling parameters are kept in the file 'config.dat' and should be self explanatory (if you do not understand the multinest parameters, the defaults are fine for most purposes). Attention should be paid to the number of live points, more live points will sample the parameter space better and so in general you should set it higher for a larger number of dimensions.


## detectors.ini
Detectors only need to be defined once in 'detectors.ini' and then you can refer to them by name in config.dat


## Running
The code can then be ran by simply typing "./DarkSearch" std out will display the progress of the sampling. 

Alternatively, using config.dat as a guide, you can create your own parameter file and pass it at runtime eg: "./DarkSearch config2.dat"

If running with mpi, you need to type "mpirun -np [number of threads] ./DarkSearch", if it complains about shared libraries try
specifying the full path to mpirun.

The results can be analysed using the mathematica package 'crediblePlot.m' in the results folder.


## Output files
Output files are automatically generated named using the 'root' parameter specified in config.dat (default is "results/DS_")

- DS_sim.dat
Contains columns of posterior probability, -2*log(likelihood) and then the sampling parameters: mass, cross section etc.

- DS_dRdE.dat
Contains columns of recoil energy and the differential event rate

The rest of the files are from multinest and are not needed, but you can find more information on them in the 
multinest README.

