/////////////////////////////
//  Likelihood functions   //
/////////////////////////////
#include <math.h>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include "directDet.h"
#include "detectors.h"
#include "detectorFunctions.h"

//natural log of Poisson dist: gives more accurate values for small probabilities (because of machine precision)
double logPoisson(double obs, double expect)
{
    if ( expect == 0. && obs == 0 )
    {
        return 0.;
    }
    else
    {
        if(obs>170)
        {
            return -expect + (double) obs * log( expect ) - ( obs * log( obs ) - obs );  // Sterling's approx.
        }	
        return -expect + (double) obs * log( expect ) - log( gsl_sf_fact( obs ) );  // or return gsl_ran_poisson_pdf (obs,expect); 
    }
}

//scales by prior and returns ranPar
double scale(double *ranPar, double x, double y, int prior)
{ 
	if( prior == 0 ) //log prior
    {
        *ranPar = pow( 10, *ranPar * ( log10(y) - log10(x) ) + log10(x) );    
    }
    else if( prior == 1 ) //linear prior
        *ranPar = *ranPar * (y - x) +  x; 
    else if( prior == 2 ) //gaussian prior
    {
        *ranPar = x + gsl_cdf_gaussian_Pinv(*ranPar, y);  
    }
    else if( prior == 3 ) //no prior
        *ranPar = x;

    return *ranPar;
}

void scaleParams(double *Cube, physicalParameters P, WIMPpars *W)
{
    
	W->Mx   = scale(&Cube[(int)P.Mx[3]],  P.Mx[0],  P.Mx[1],  (int)P.Mx[2]);	
    for(int i=1;i<16;i++)
        W->coeff[i] = scale(&Cube[(int)P.coeff[i][3]], P.coeff[i][0], P.coeff[i][1], (int)P.coeff[i][2]);
       
    W->rho  = scale(&Cube[(int)P.rho[3]], P.rho[0], P.rho[1], (int)P.rho[2]);
    W->v0   = scale(&Cube[(int)P.v0[3]],  P.v0[0],  P.v0[1],  (int)P.v0[2]);
    W->vesc = scale(&Cube[(int)P.vesc[3]],P.vesc[0],P.vesc[1],(int)P.vesc[2]);
    W->vSp  = scale(&Cube[(int)P.vSp[3]],P.vSp[0],P.vSp[1],(int)P.vSp[2]);
    W->vEp  = scale(&Cube[(int)P.vEp[3]],P.vEp[0],P.vEp[1],(int)P.vEp[2]);
    
    W->spin = P.spin;
    
}

//likelihood function for binned data
double logLikelihood( WIMPpars W, detector *dets, int ndets, physicalParameters P, int b)
{

	//Calculate log-likelihood
	double predicted, background, l;
	double loglike = 0;
    double Er_min, Er_max;

    //loop over detectors
    for(int j=0; j<ndets; j++)
    {   
        
        //loop over recoil energy bins
        for(int i=0; i<dets[j].nbins; i++)			
        {
            Er_min = (double)i*dets[j].binW + dets[j].ErL;
            Er_max = (double)(i+1)*dets[j].binW + dets[j].ErL;
            
            background = b * BgRate(dets[j], Er_min, Er_max) * dets[j].exposure;
            predicted = intRate( Er_min, Er_max, W, dets[j], P) * dets[j].exposure; 
            
            l = logPoisson( dets[j].binnedData[i], predicted+background);
            loglike += l;
            
        } 
        
    }

	return loglike;
}

//binless likelihood function
double logLikelihoodBinless( WIMPpars W, detector *dets, int ndets, physicalParameters P, int b)
{

	//Calculate log-likelihood
	double predicted, background;
	double loglike = 0;
    
    //loop over detectors
    for(int j=0; j<ndets; j++)
    {   
    
        //total expected background
        background = b * BgRate(dets[j], dets[j].ErL, dets[j].ErU) * dets[j].exposure;
        predicted = intRate( dets[j].ErL, dets[j].ErU, W, dets[j], P) * dets[j].exposure; 
        loglike += logPoisson( dets[j].nEvents, predicted+background);
        
        //loop over events
        for (int i=0; i<dets[j].nEvents; i++)
        {
            //unfinshed likelihood
            loglike += dets[j].exposure/predicted * diffRate( dets[j].unbinnedData[i], W, dets[j], P);
        }
            
    }

	return loglike;
}

//log likelihood function for differential rate dN from WIMPs (background not included)
//This is the one that gets passed to MultiNest
void LogLikedN(double *Cube, int &ndim, int &npars, double &lnew, long &pointer)    
{   

	//get pointer in from MultiNest 
    parameterList pL = *(parameterList *) pointer;
    
    //WIMP pars for this point in the parameter space
    WIMPpars Wcube;
	scaleParams( Cube, pL.p, &Wcube);
	
    if(pL.binlessL==1)
    {
        lnew = logLikelihoodBinless( Wcube, pL.detectors, pL.ndet, pL.p, 1);
    }
    else
    {
        lnew = logLikelihood( Wcube, pL.detectors, pL.ndet, pL.p, 1);
    }
}

