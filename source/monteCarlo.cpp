#include <cmath>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include "parameterStruct.h"
#include "directDet.h"
#include "detectors.h"
#include "detectorFunctions.h"

 extern const gsl_rng_type * T;
 extern gsl_rng * r; 

int SEED = 0;

//Generates count number of random events, distributed according to dN/dE for parameters in cube, for detector det (0 or 1, corresponding to order in sampling.par), recoil energies [keV] are stored in MCdata
void generateUnbinnedData(WIMPpars W, physicalParameters P, detector *det, int b)
{

    double predicted = intRate( det->ErL, det->ErU, W, *det, P) * det->exposure; 
    double background= BgRate(*det, det->ErL, det->ErU) * det->exposure; 
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, (int)time(NULL)+ SEED++);
    
    int count = predicted+background;
    
    det->unbinnedData = new double[count];
    
    double norm = diffRate( det->ErL, W, *det, P);
    int i = 0;
    double Er,y;
    while ( i < count)
    {
        Er = gsl_rng_uniform(r)*(det->ErU-det->ErL) + det->ErL;
        y = gsl_rng_uniform(r);
        if( diffRate( Er, W, *det, P)/norm  > y )
        {
            det->unbinnedData[i] = Er;
            i++;
        }
    }
}

void generateBinnedData(WIMPpars W, physicalParameters P, detector *det, int b)
{
    double predicted, background; 
    double Er_min, Er_max;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, (int)time(NULL)+ SEED++);
    
    predicted = (intRate( det->ErL, det->ErU, W, *det, P) + BgRate(*det, det->ErL, det->ErU)) * det->exposure; 

    det->nbins = floor( sqrt(predicted)/2 );  
    if(det->nbins == 0)
        det->nbins = 1;
    det->binW = ( det->ErU - det->ErL ) / ( (double) det->nbins);
    
    det->binnedData = new double[det->nbins];
    for(int i=0; i<det->nbins; i++)
    {
        Er_min = (double)i*det->binW+det->ErL;
        Er_max = (double)(i+1)*det->binW+det->ErL;
         
        background = b * BgRate(*det, Er_min, Er_max) * det->exposure;
        predicted = intRate( Er_min, Er_max, W, *det, P) * det->exposure; 
        
        if( W.asimov == 1) 
            det->binnedData[i] = gsl_ran_poisson(r,predicted+background);
        else
            det->binnedData[i] = predicted + background;            
        
        det->nEvents += det->binnedData[i];
    }
     
}

