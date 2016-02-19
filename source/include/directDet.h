#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double intWIMPrate(double Er_min, double Er_max, WIMPpars *W, detector *det);

double diffWIMPrate(double Er, WIMPpars *W, detector *det);

double smearedDiffWIMPrate(double Er, WIMPpars *W, detector *det);
