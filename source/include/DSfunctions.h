//Header for Darth Sidious direct detection related functions

//C includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <assert.h>

//GNU Science lib includes
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>

//Other includes
#include "DEIntegrator.h"
#include "multinest.h"

//Darth Sidious includes
#include "directDet.h"
#include "likelihood.h"
#include "monteCarlo.h"
#include "parameterStruct.h"

using namespace std;

#ifdef MPI
    #include "mpi.h"
#endif

//needed by multinest for monitoring progress, I don't use it
void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void **pointer) {}

int getSamplingPars(parameterList *pL, char *filename);

int twoDbin(char *filename, int nbin);           //2D binning of the output from multinest (marginalisation)

int bin(double *data, int lenData, double min, double max); //return number of events between min and max