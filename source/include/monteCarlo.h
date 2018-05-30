#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int generateBinnedData(WIMPpars *W, detector *det, int nbins, int b, int simSeed);

int generateUnbinnedData(WIMPpars *W, detector *det, int b, int simSeed);

int generateTimeBinnedData(WIMPpars *W, detector *det, int nbins, int b, int simSeed);
