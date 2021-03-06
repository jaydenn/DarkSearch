//Header for Darth Sidious related functions

//C includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <time.h>

//GNU Science lib includes
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>

//Other includes
#include "DEIntegrator.h"
#include "multinest.h"

//Dark Search includes
#include "directDet.h"
#include "detectors.h"
#include "likelihood.h"
#include "monteCarlo.h"
#include "nuBackground.h"
#include "fileIO.h"

#ifdef MPI
    #include "mpi.h"
#endif

//set up random number generator for Monte Carlo simulation
  const gsl_rng_type * T;
  gsl_rng * r;

//needed by multinest for monitoring progress, I don't use it
//void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void **pointer) {}
void dumper(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void **) {}

