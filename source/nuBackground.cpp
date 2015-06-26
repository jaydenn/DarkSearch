#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "formFactorSI.h"
#include "DEIntegrator.h"
#include "detectorFunctions.h"
#include "nuFluxData.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>


struct intParamStruct {

    double Er;
    double A;
    gsl_interp *flux;
    gsl_interp_accel *accel;
    
};

double integrand(double Enu, void *pars)
{
    intParamStruct *intPars = (intParamStruct*)pars;
    return gsl_interp_eval(intPars->flux,flux_E,flux_N,Enu,intPars->accel) * (1 - 938000.*intPars->A*intPars->Er/(2*Enu*Enu)) ;
}

//returns rate of neutrino events in detector D per (tonne.year.keV) at the specified Er
double dNnudE(detector det, double Er)  						  
{

    gsl_interp *flux = gsl_interp_alloc(gsl_interp_linear,flux_points);
    gsl_interp_init(flux,flux_E,flux_N,flux_points);
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    
    double EnuMax = flux_E[flux_points-1];
    double EnuMin = sqrt( Er * 938000.*det.AM / 2);
    if (EnuMin < 100)
    {
        std::cout << "EnuMin = " << Er << " " << EnuMin << " <100keV, bg rate may be inaccurate" << std::endl;
        EnuMin=100;
    }
    
    intParamStruct intPars; //intParamStruct();
    intPars.A = det.AM;
    intPars.Er = Er;
    intPars.flux = flux;
    intPars.accel = accel;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    F.function = &integrand;
    F.params = &intPars;
    double integral,absErr;
    
    int key = 3;
    int limit = 1000;
    
    gsl_integration_qag(&F, EnuMin, EnuMax, 10, 1e-4, limit, key, w, &integral, &absErr);
    
    gsl_interp_free(flux);
    gsl_integration_workspace_free(w);
    
    return 7.484e-8 * pow( ffactorSI( det.AM, Er),2) * pow(det.AM-1.04*det.isoZ[0],2) * integral;
}

