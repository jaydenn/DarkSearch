#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include "parameterStruct.h"
#include "directDet.h"
#include "detectors.h"
#include "detectorFunctions.h"

int SEED=0;

//Generates random number of randomly distributed according to dN/dE for parameters in cube, for detector det (0 or 1, corresponding to order in sampling.par), recoil energies [keV] are stored in MCdata
int generateUnbinnedData(WIMPpars *W, detector *det, int b, int simSeed)
{

    double signal = intWIMPrate( det->ErL, det->ErU, W, det) * det->exposure;   
    double background= b*intBgRate(*det, det->ErL, det->ErU) * det->exposure; 
     
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, simSeed + SEED++);
    
    det->nEvents = gsl_ran_poisson(r,signal+background);
    
    det->unbinnedData = new double[(int)det->nEvents];
    
    double norm = diffWIMPrate( det->ErL, W, det) + b*diffBgRate(*det,det->ErL);
    int i = 0;
    double Er,y;
    while ( i < det->nEvents)
    {
        Er = gsl_rng_uniform(r)*(det->ErU-det->ErL) + det->ErL;
        y = gsl_rng_uniform(r);
        if( (diffWIMPrate( det->ErL, W, det) + b*diffBgRate(*det,det->ErL))/norm  > y )
        {
            det->unbinnedData[i] = Er;
            i++;
        }
    }
    
    gsl_rng_free(r);
    return 0;
}

int generateBinnedData(WIMPpars *W, detector *det, int nbins, int b, int simSeed)
{

    double Er_min, Er_max;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, simSeed + SEED++);
    
    //total signal and background, for setting bin size
    double signal = intWIMPrate( det->ErL, det->ErU, W, det) * det->exposure;
    double background = b * intBgRate(*det, det->ErL, det->ErU) * det->exposure; 

    det->nbins = nbins;

    det->binW = ( det->ErU - det->ErL ) / ( (double) det->nbins);
    
    try
    {
        det->binnedData = new double[det->nbins];
    }
    catch (std::bad_alloc& ba)
    {
      std::cerr << "bad_alloc caught: " << ba.what() << std::endl << "you requested: " << det->nbins << " doubles" <<std::endl;
      return 1;
    }
    
    char filename[100];
    std::ofstream datfile;
    FILE* data;
    char *ret;
    char temp[100];
    
    sprintf(filename,"%s_events.dat",det->name);
    if ( W->asimov == 2 )
        data = fopen(filename,"r");
    else
        datfile.open(filename,ios::out);

    for(int i=0; i<det->nbins; i++)
    {
        Er_min = (double)i*det->binW+det->ErL;
        Er_max = (double)(i+1)*det->binW+det->ErL;
        
        background = b * intBgRate(*det, Er_min, Er_max) * det->exposure;
        signal = intWIMPrate( Er_min, Er_max, W, det) * det->exposure; 
        
        if ( W->asimov == 2 )
        {
            ret = fgets(temp,100,data);
            sscanf(temp,"%*lf %*lf %lf",&(det->binnedData[i]));
        }
        else if( W->asimov == 1)
        {
            det->binnedData[i] = gsl_ran_poisson(r,signal+background);
            datfile << Er_min << " " << Er_max << " " << det->binnedData[i] << std::endl;
        }
        else
        {
            det->binnedData[i] = signal + background;            
            datfile << Er_min << " " << Er_max << " " << det->binnedData[i] << std::endl;
        }
        
        det->nEvents += det->binnedData[i];
    }

    datfile.close();
    gsl_rng_free(r);
    return 0;
     
}

int generateTimeBinnedData(WIMPpars *W, detector *det, int nbins, int b, int simSeed)
{

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, simSeed + SEED++);
    
    //total signal and background, for setting bin size
    double signal;// = intWIMPrate( det->ErL, det->ErU, W, det) * det->exposure;
    double background;// = b * intBgRate(*det, det->ErL, det->ErU) * det->exposure; 

    det->nbins = nbins;

    det->binW = 1 / ( (double) det->nbins);
    
    try
    {
        det->binnedData = new double[det->nbins];
    }
    catch (std::bad_alloc& ba)
    {
      std::cerr << "bad_alloc caught: " << ba.what() << std::endl << "you requested: " << det->nbins << " doubles" <<std::endl;
      return 1;
    }
    
    char filename[100];
    std::ofstream datfile;
    FILE* data;
    char *ret;
    char temp[100];
    
    sprintf(filename,"%s_events.dat",det->name);
    if ( W->asimov == 2 )
        data = fopen(filename,"r");
    else
        datfile.open(filename,ios::out);
    
    double T_min,T_max;
    for(int i=0; i<det->nbins; i++)
    {
        T_min = (double)i*det->binW;
        T_max = (double)(i+1)*det->binW;
        
        //background = b * intBgRate(*det, Er_min, Er_max) * det->exposure;
        signal = intWIMPrateT( det->ErL, det->ErU, T_min, T_max, W, det) * det->exposure; 
        
        if ( W->asimov == 2 )
        {
            ret = fgets(temp,100,data);
            sscanf(temp,"%*lf %*lf %lf",&(det->binnedData[i]));
        }
        else if( W->asimov == 1)
        {
            det->binnedData[i] = gsl_ran_poisson(r,signal+background);
            datfile << T_min << " " << T_max << " " << det->binnedData[i] << std::endl;
        }
        else
        {
            det->binnedData[i] = signal + background;            
            datfile << T_min << " " << T_max << " " << det->binnedData[i] << std::endl;
        }
        std::cout << i << " " << T_min << " " << det->binnedData[i] << std::endl;
        det->nEvents += det->binnedData[i];
    }

    datfile.close();
    gsl_rng_free(r);
    return 0;
     
}

