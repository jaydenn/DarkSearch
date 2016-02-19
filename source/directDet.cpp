#include <cmath>
#include <iostream>
#include "DEIntegrator.h"
#include "parameterStruct.h"
#include "detectorFunctions.h"
#include "detectors.h"
#include "isoRates.h"

//returns direct detection rate /t/year/keV evaluated at recoil Er without smearing for detector resolution
double unsmearedDiffWIMPrate(double Er, WIMPpars *W, detector *D)						 
{

    //contributions for each isotope
    double totRate = 0;
    for(int i=0; i<D->nIso; i++)
    {
        totRate += D->isoFrac[i] * isoRate(Er, W, D->isoA[i], D->isoZ[i]);              
    }

    return totRate * detEff(Er,D->eff);

}

class convolveIntegral
{
public:
    double E;
    WIMPpars *W;
    detector *D;
    double operator()(double Er) const
    {
        return unsmearedDiffWIMPrate( Er, W, D) * 1/(sqrt(2*M_PI)*detRes(E,D->res)) * exp( -pow(E-Er,2)/(2*pow(detRes(E,D->res),2)) );
    }
};

double smearedDiffWIMPrate(double Er, WIMPpars *W, detector *D) 
{
    convolveIntegral conInt;   //create instance for recoil energy integration
    conInt.W = W;              //set the integral parameters
    conInt.D = D;
    conInt.E = Er;
    return DEIntegrator<convolveIntegral>::Integrate(conInt,0,100,1e-5);
}

//returns direct detection rate /t/year/keV evaluated at recoil Er 
double diffWIMPrate(double Er, WIMPpars *W, detector *D)						 
{
    if(detRes(Er,D->res) > 1e-10)                    //To smear, or not to smear?
        return smearedDiffWIMPrate(Er, W, D);
    else
        return unsmearedDiffWIMPrate(Er, W, D);
}

class ErIntegral
{
public:
    WIMPpars *W;
    detector *D;
    double operator()(double Er) const
    {
        return diffWIMPrate( Er, W, D);
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

