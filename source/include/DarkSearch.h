//Header for Darth Sidious related functions

//C includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <assert.h>
#include <time.h>

//GNU Science lib includes
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>

//Other includes
#include "DEIntegrator.h"
#include "multinest.h"

//Darth Sidious includes
#include "directDet.h"
#include "detectors.h"
#include "likelihood.h"
#include "monteCarlo.h"
#include "nuBackground.h"

#ifdef MPI
    #include "mpi.h"
#endif

using namespace std;

//set up random number generator for Monte Carlo simulation
  const gsl_rng_type * T;
  gsl_rng * r;

//needed by multinest for monitoring progress, I don't use it
//void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void **pointer) {}
void dumper(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void **) {}

int twoDbin(char *filename, int x, int y, int nbin);           //2D binning of the output from multinest (marginalisation)

int bin(double *data, int lenData, double min, double max); //return number of events between min and max

int getSamplingPars(parameterList *pL, char *filename);

