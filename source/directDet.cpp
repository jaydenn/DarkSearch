#include <cmath>
#include <iostream>
#include "DEIntegrator.h"
#include "parameterStruct.h"
#include "detectorFunctions.h"
#include "detectors.h"
#include "isoRates.h"

//returns direct detection rate /t/year/keV evaluated at recoil Er
double diffRate(double Er, WIMPpars W, detector D, physicalParameters P)						 
{

    //contributions for each isotope
    double totRate = 0;
    for(int i=0; i<D.nIso; i++)
    {
        totRate += D.isoFrac[i] * isoRate(Er, W, P, D.isoA[i], D.isoZ[i]);              
    }
   
    return totRate * detEff(Er,D.eff);

}

class convolveIntegral
{
public:
    double E;
    WIMPpars W;
    detector D;
    physicalParameters P;
    double operator()(double Er) const
    {
        return diffRate( Er, W, D, P) * 1/(sqrt(2*M_PI)*detRes(E,D.res)) * exp( -pow(E-Er,2)/(2*pow(detRes(E,D.res),2)) );
    }
};

double smearedDiffRate(double Er, WIMPpars W, detector D, physicalParameters P) 
{
    convolveIntegral conInt;   //create instance for recoil energy integration
    conInt.W = W;              //set the integral parameters
    conInt.D = D;
    conInt.P = P;
    conInt.E = Er;
    return DEIntegrator<convolveIntegral>::Integrate(conInt,0,100,1e-5);
  
}


class smearedErIntegral
{
public:
    WIMPpars W;
    detector D;
    physicalParameters P;
    double operator()(double Er) const
    {
        return smearedDiffRate( Er, W, D, P);
    }
};

class ErIntegral
{
public:
    WIMPpars W;
    detector D;
    physicalParameters P;
    double operator()(double Er) const
    {
        return diffRate( Er, W, D, P);
    }
};

//returns total integrated event rate per tonne/year in detector with recoil Er_min < Er < Er_max
double intRate(double Er_min, double Er_max, WIMPpars W, detector D, physicalParameters P)						  
{
    
    /*if(D.res > 10)                    //To smear, or not to smear?
    {
        smearedErIntegral sErInt;      //create instance for smeared recoil energy integration
        sErInt.W = W;
        sErInt.D = D;
        sErInt.P = P;

        return DEIntegrator<smearedErIntegral>::Integrate(sErInt,Er_min,Er_max,1e-4);
    }
    else
    {*/
        ErIntegral erInt; 		       //create instance for recoil energy integration
        erInt.W = W;
        erInt.D = D;
        erInt.P = P;
        
        return DEIntegrator<ErIntegral>::Integrate(erInt,Er_min,Er_max,1e-4);      
    //}
}
