#include <cmath>
#include <iostream>
#include "DEIntegrator.h"
#include "parameterStruct.h"
#include "detectorFunctions.h"
#include "detectors.h"
#include "isoRates.h"

//returns direct detection rate /t/year/keV evaluated at recoil Er without smearing for detector resolution
double unsmearedDiffWIMPrate(double Er, WIMPpars *W, detector *D, double T)						 
{

    //contributions for each isotope
    double totRate = 0;
    for(int i=0; i<D->nIso; i++)
    {
        totRate += D->isoFrac[i] * isoRateT(Er, W, D->isoA[i], D->isoZ[i], T);              
    }

    return totRate * detEff(Er,D->eff);

}

class convolveIntegral
{
public:
    double E;
    WIMPpars *W;
    detector *D;
    double T;
    double operator()(double Er) const
    {
        return unsmearedDiffWIMPrate( Er, W, D, T) * 1/(sqrt(2*M_PI)*detRes(E,D->res)) * exp( -pow(E-Er,2)/(2*pow(detRes(E,D->res),2)) );
    }
};

double smearedDiffWIMPrate(double Er, WIMPpars *W, detector *D, double T) 
{
    convolveIntegral conInt;   //create instance for recoil energy integration
    conInt.W = W;              //set the integral parameters
    conInt.D = D;
    conInt.E = Er;
    conInt.T = T;
    double ErLow = Er-3*detRes(Er,D->res);
    if(ErLow < 0)
        ErLow = 0;
    return DEIntegrator<convolveIntegral>::Integrate(conInt,ErLow,Er+3*detRes(Er,D->res),1e-4);
}

//returns direct detection rate /t/year/keV evaluated at recoil Er 
double diffWIMPrateT(double Er, WIMPpars *W, detector *D, double T)						 
{
    if(detRes(Er,D->res) > 1e-10)                    //To smear, or not to smear?
        return smearedDiffWIMPrate(Er, W, D, T);
    else
        return unsmearedDiffWIMPrate(Er, W, D, T);
}

double diffWIMPrate(double Er, WIMPpars *W, detector *D)						 
{
    return diffWIMPrateT(Er,W,D,0);
}

class ErIntegral
{
public:
    WIMPpars *W;
    detector *D;
    double T;
    double operator()(double Er) const
    {
        return diffWIMPrateT( Er, W, D, T);
    }
};

//returns total integrated event rate per tonne/year in detector with recoil Er_min < Er < Er_max
double intWIMPrate(double Er_min, double Er_max, WIMPpars *W, detector *D)						  
{

    ErIntegral erInt; 		       //create instance for recoil energy integration
    erInt.W = W;
    erInt.D = D;
    
    return DEIntegrator<ErIntegral>::Integrate(erInt,Er_min,Er_max,1e-4);
    
}

class TimeIntegral
{
public:
    WIMPpars *W;
    detector *D;
    ErIntegral *erInt;
    double Er_min, Er_max;
    double operator()(double T) const
    {
        erInt->T = T;
        return DEIntegrator<ErIntegral>::Integrate(*erInt,Er_min,Er_max,1e-4);
    }
};


//returns total integrated event rate per tonne/year in detector with recoil Er_min < Er < Er_max and T_min < T < T_max
double intWIMPrateT(double Er_min, double Er_max, double T_min, double T_max, WIMPpars *W, detector *D)						  
{

    ErIntegral erInt; 		       //create instance for recoil energy integration
    erInt.W = W;
    erInt.D = D;
    TimeIntegral tInt; 		       //create instance for recoil energy integration
    tInt.W = W;
    tInt.D = D;
    tInt.erInt = &erInt;
    tInt.Er_min = Er_min;
    tInt.Er_max = Er_max;
    
    return DEIntegrator<TimeIntegral>::Integrate(tInt,T_min,T_max,1e-3);
    
}

