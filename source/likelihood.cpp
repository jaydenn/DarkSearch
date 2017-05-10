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
    if ( expect > 0. && obs > 0. )
        return -expect + obs * log( expect ) - gsl_sf_lngamma( obs+1 );
    else
        return -1E299;
}

//scales by prior and returns ranPar
double scale(double *ranPar, double x, double y, int prior)
{ 
	if( prior == 0 ) //log prior
        *ranPar = pow( 10, *ranPar * ( log10(y) - log10(x) ) + log10(x) );    
    else if( prior == 1 ) //linear prior
        *ranPar = *ranPar * (y - x) +  x; 
    else if( prior == 2 ) //gaussian prior
        *ranPar = x + gsl_cdf_gaussian_Pinv(*ranPar, y);  
    else if( prior == 3 ) //no prior
        return x;

    return *ranPar;
}

void scaleParams(double *Cube, reconstructionParameters P, WIMPpars *W)
{
    
	W->Mx   = scale(&Cube[(int)P.Mx[3]],  P.Mx[0],  P.Mx[1],  (int)P.Mx[2]);

    if( ((int) P.spin[2])==3)
        W->spin = P.spin[0];    
    else
    {
        if ( Cube[(int)P.spin[3]] < 0.6666 )
        {
            if ( Cube[(int)P.spin[3]] < 0.3333 )
                W->spin = Cube[(int)P.spin[3]] = 0;
            else
                W->spin = Cube[(int)P.spin[3]] = 0.5;
        }
        else
            W->spin = Cube[(int)P.spin[3]] = 1;
    }
    
    for(int i=1;i<16;i++)
    {
        if(P.isv == 0) //proton and neutron values fixed to eachother
        {
            W->coeffn[i] = scale(&Cube[(int)P.coeffn[i][3]], P.coeffn[i][0], P.coeffn[i][1], (int)P.coeffn[i][2]);
            W->coeffp[i] = W->coeffn[i];
        }
        else if (P.isv == 1) //proton and neutron values allowed to be different
        {
            W->coeffn[i] = scale(&Cube[(int)P.coeffn[i][3]], P.coeffn[i][0], P.coeffn[i][1], (int)P.coeffn[i][2]);
            W->coeffp[i] = scale(&Cube[(int)P.coeffp[i][3]], P.coeffp[i][0], P.coeffp[i][1], (int)P.coeffp[i][2]);
        }
        else if (P.isv == 2) //proton scanned, neutron set to zero
        {
            W->coeffn[i] = 0;
            W->coeffp[i] = scale(&Cube[(int)P.coeffp[i][3]], P.coeffp[i][0], P.coeffp[i][1], (int)P.coeffp[i][2]);
        }
        else if (P.isv == 3) //neutron scanned, proton set to zero
        {
            W->coeffn[i] = scale(&Cube[(int)P.coeffn[i][3]], P.coeffn[i][0], P.coeffn[i][1], (int)P.coeffn[i][2]);
            W->coeffp[i] = 0;
        }
    }
       
    W->rho  = scale(&Cube[(int)P.rho[3]], P.rho[0], P.rho[1], (int)P.rho[2]);

    W->vDindex=P.vDindex;
    
    if(P.vDindex>1) //if using polynomials for velocity dist
    {
        for(int i=1;i<=P.vDindex;i++)
            W->vLa[i] = scale(&Cube[(int)P.vLa[i][3]],  P.vLa[i][0],  P.vLa[i][1],  (int)P.vLa[i][2]);
    }
    else
    {
        W->v0   = scale(&Cube[(int)P.v0[3]], P.v0[0],  P.v0[1],  (int)P.v0[2]);
        W->vesc = scale(&Cube[(int)P.vesc[3]], P.vesc[0],P.vesc[1],(int)P.vesc[2]);
        W->vSp  = scale(&Cube[(int)P.vSp[3]], P.vSp[0],P.vSp[1],(int)P.vSp[2]);
        W->vEp  = scale(&Cube[(int)P.vEp[3]], P.vEp[0],P.vEp[1],(int)P.vEp[2]);
    }
}

//likelihood function for binned data
double logLikelihood( WIMPpars *W, detector *dets, int ndets, int b)
{

	//Calculate log-likelihood
	double signal, background, l;
	double loglike = 0;
    double Er_min, Er_max;

    //loop over detectors
    for(int j=0; j<ndets; j++)
    {   
        
        //loop over recoil energy bins
        for(int i=0; i<dets[j].nbins; i++)			
        {
            //set bin limits
            Er_min = (double)i*dets[j].binW + dets[j].ErL;
            Er_max = (double)(i+1)*dets[j].binW + dets[j].ErL;
            
            background = b * intBgRate(dets[j], Er_min, Er_max) * dets[j].exposure;
            signal = intWIMPrate(Er_min, Er_max, W, &(dets[j])) * dets[j].exposure; 

            l = logPoisson( dets[j].binnedData[i], signal+background+1e-99);
            loglike += l;
            
        } 
        
    }

	return loglike;
}

//binless likelihood function
double logLikelihoodBinless( WIMPpars *W, detector *dets, int ndets, int b)
{

	//Calculate log-likelihood
	double signal, background;
	double loglike = 0;
    
    //loop over detectors
    for(int j=0; j<ndets; j++)
    {   
    
        //total expected events
        background = b * intBgRate( dets[j], dets[j].ErL, dets[j].ErU) * dets[j].exposure;
        signal  =   intWIMPrate(dets[j].ErL, dets[j].ErU, W, &(dets[j])) * dets[j].exposure; 
        loglike += logPoisson( dets[j].nEvents, signal+background);
            
        //loop over events
        for (int i=0; i<dets[j].nEvents; i++)
            loglike += log( dets[j].exposure * ( diffWIMPrate( dets[j].unbinnedData[i], W, &(dets[j])) + diffBgRate(dets[j],dets[j].unbinnedData[i]) ) / (signal + background) );
            
    }

	return loglike;
}

//log likelihood function for differential rate dN from WIMPs (background not included)
//This is the one that gets passed to MultiNest
void LogLikedN(double *Cube, int &ndim, int &npars, double &lnew, long &pointer)    
{   

    //get pointer in from MultiNest 
    parameterList *pL = (parameterList *) pointer;
    
    //WIMP pars for this point in the parameter space
    WIMPpars Wcube;
	scaleParams( Cube, pL->p, &Wcube);
	
    if(pL->binlessL==1)
    {
        lnew = logLikelihoodBinless( &Wcube, pL->detectors, pL->ndet, 1);
    }
    else
    {
        lnew = logLikelihood( &Wcube, pL->detectors, pL->ndet, 1);
    }
    
}

