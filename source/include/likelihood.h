#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double logPoisson(int obs, double expect);

void scaleParams(double *Cube, reconstructionParameters P, WIMPpars W);

//likelihood function for binned data
double logLikelihood( WIMPpars *W, detector *dets, int ndets, reconstructionParameters P, int b);

//likelihood function for unbinned data
double logLikelihoodBinless( WIMPpars *W, detector *dets, int ndets, reconstructionParameters P, int b);

//log likelihood function for passing the above to Multinest
void LogLikedN(double *Cube, int &ndim, int &npars, double &lnew, long &pointer);

