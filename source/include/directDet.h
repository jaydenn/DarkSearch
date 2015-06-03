#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double intRate(double Er_min, double Er_max, WIMPpars W, detector det, physicalParameters P);

double diffRate(double Er, WIMPpars W, detector det, physicalParameters P);

double smearedDiffRate(double Er, WIMPpars W, detector det, physicalParameters P);
