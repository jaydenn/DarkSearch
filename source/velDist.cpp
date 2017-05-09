#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_legendre.h>
#include "DEIntegrator.h"

double fLegendre(double v, double *a)
{
    double Ps[5]; double sum=0;
    gsl_sf_legendre_Pl_array(4, 2*v/1000.0 - 1, Ps);
  
    for( int i=0; i < 5; i++)
        sum -= a[i]*Ps[i];

    return v*v*exp( sum );
}

class velocityIntegralNorm
{
public:
    double *a;
    double operator()(double v) const
    {
        return fLegendre( v, a);
    }
};

class velocityIntegralonV
{
public:
    double *a;
    double operator()(double v) const
    {
        return fLegendre( v, a)/v;
    }
};

class velocityIntegralvSq
{
public:
    double *a;
    double operator()(double v) const
    {
        return fLegendre( v, a)*v;
    }
};

double G(double v0, double ve, double vesc, double vmin, int index, double *a)
{
    
    switch (index)
    {
        case 0:
        {
            return ( erfc( (ve - vmin) / v0 ) + erfc( (ve + vmin) / v0 ) ) / (2*ve);
        }
         
        case 1:
        {
            double x = 0;
            
            if(-(ve/v0) + vesc/v0 - vmin/v0 > 0) 
                x += v0*((-4*ve*(-3 + pow(ve,2)/pow(v0,2) - (3*pow(vesc,2))/pow(v0,2) + (3*pow(vmin,2))/pow(v0,2)))/v0 + 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*(-erf(ve/v0 - vmin/v0) - erf(ve/v0 + vmin/v0)))/    (2.*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0))); 
            if(ve/v0 + vesc/v0 - vmin/v0 > 0 && ve/v0 - vesc/v0 + vmin/v0 > 0)
                x += (v0*(2*(ve/v0 + vesc/v0 - vmin/v0)*(3 + (ve/v0 + vesc/v0 - vmin/v0)*(-(ve/v0) + (2*vesc)/v0 + vmin/v0)) + 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*(-erf(vesc/v0) - erf(ve/v0 - vmin/v0))))/(2.*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0)));
        
            return x/v0;
        }
        
        case 2:
        {
            velocityIntegralonV vInt;      //create instance for smeared recoil energy integration
            velocityIntegralNorm vIntN;
            vInt.a = a;
            vIntN.a = a;
            
            double norm = DEIntegrator<velocityIntegralNorm>::Integrate(vIntN,0,1000,1e-5);
            return DEIntegrator<velocityIntegralonV>::Integrate(vInt,vmin,1000,1e-5)/norm;
        }
        
        default:
        {
            std::cout << "vDist not found: " << index << std::endl;
            return(0);
        }
    }
}

double Gsq(double v0, double ve, double vesc, double vmin, int index, double *a)
{
    switch (index)
    {
        case 0:
        {
            return ( erfc( (ve - vmin) / v0 ) + erfc( (ve + vmin) / v0 ) ) / (2*ve);
        }
        
        case 1:
        {
             double x = 0;
            if(-(ve/v0) + vesc/v0 - vmin/v0 > 0) 
                x += -(v0*(-30*exp(pow(vesc,2)/pow(v0,2) + pow(-(ve/v0) + vmin/v0,2))*(-(ve/v0) + vmin/v0) + 30*exp(pow(vesc,2)/pow(v0,2) + pow(ve/v0 + vmin/v0,2))*(ve/v0 + vmin/v0) +  (4*exp(2*(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2)))*ve*(pow(ve,4)/pow(v0,4) - (10*pow(ve,2)*(1 + pow(vesc,2)/pow(v0,2)))/pow(v0,2) - 15*(2 + (2*pow(vesc,2))/pow(v0,2) + pow(vesc,4)/pow(v0,4) - pow(vmin,4)/pow(v0,4))))/v0 + 15*exp(pow(vesc,2)/pow(v0,2) + 2*(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2)))*sqrt(M_PI)*(1 + (2*pow(ve,2))/pow(v0,2))*(erf(ve/v0 - vmin/v0) + erf(ve/v0 + vmin/v0))))/(20.*exp(2*(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2)))*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0)));
            if(ve/v0 + vesc/v0 - vmin/v0 > 0 && ve/v0 - vesc/v0 + vmin/v0 > 0)
                x +=    (exp(-(pow(ve,2)/pow(v0,2)) - pow(vmin,2)/pow(v0,2))*v0*(2*exp(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2))*(-(pow(ve,5)/pow(v0,5)) + (15*vesc)/v0 + (10*pow(vesc,3))/pow(v0,3) + (4*pow(vesc,5))/pow(v0,5) + (10*pow(ve,3)*(1 + pow(vesc,2)/pow(v0,2)))/pow(v0,3) - (10*pow(vmin,3))/pow(v0,3) - (10*pow(vesc,2)*pow(vmin,3))/pow(v0,5) + (6*pow(vmin,5))/pow(v0,5) + (10*pow(ve,2)*((3*vesc)/v0 + (2*pow(vesc,3))/pow(v0,3) + pow(vmin,3)/pow(v0,3)))/pow(v0,2) + (15*ve*(2 + (2*pow(vesc,2))/pow(v0,2) + pow(vesc,4)/pow(v0,4) - pow(vmin,4)/pow(v0,4)))/v0) + 15*exp(pow(vesc,2)/pow(v0,2))*(-2*exp((2*ve*vmin)/pow(v0,2))*(ve/v0 + vmin/v0) + exp(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2))*sqrt(M_PI)*(1 + (2*pow(ve,2))/pow(v0,2))*(-erf(vesc/v0) - erf(ve/v0 - vmin/v0)))))/(20.*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0)));
        
            return x*v0;
        }
        
        case 2:
        {
            velocityIntegralvSq vInt;      //create instance for smeared recoil energy integration
            velocityIntegralNorm vIntN;
            vInt.a = a;
            vIntN.a = a;
            
            double norm = DEIntegrator<velocityIntegralNorm>::Integrate(vIntN,0,1000,1e-5);
            return DEIntegrator<velocityIntegralvSq>::Integrate(vInt,vmin,1000,1e-5)/norm;
        }
        default:
        {
            std::cout << "vDist not found" << std::endl;
            return(0);
        }
    }
}
