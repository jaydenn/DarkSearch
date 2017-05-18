#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_legendre.h>
#include "DEIntegrator.h"
#include <gsl/gsl_spline.h>

#define VMAX 0.0033333

gsl_spline *FonVspline = gsl_spline_alloc(gsl_interp_linear,1000);
gsl_spline *FxVspline = gsl_spline_alloc(gsl_interp_linear,1000);
gsl_interp_accel *FonVaccel = gsl_interp_accel_alloc();
gsl_interp_accel *FxVaccel = gsl_interp_accel_alloc();

void chebyshev_array(int i, double x, double *result)
{
    result[0] = 1;
    result[1] = x;
    result[2] = 2*x*x - 1;
    if(i > 2)
    {
        result[3] = 4*x*x*x - 3*x;
        if(i > 3)
        {
            result[4] = 1-8*pow(x,2) + 8*pow(x,4);
            if(i > 4)
            {
                result[5] = 5*x - 20*pow(x,3) + 16*pow(x,5);
                if(i > 5)
                {
                    result[6] = -1 + 18*pow(x,2) - 48*pow(x,4) + 32*pow(x,6);
                    if(i > 6)
                        result[7] = -7*x + 56*pow(x,3) - 112*pow(x,5) + 64*pow(x,7);
                }
            }
        }
    }                   
}

// only normalized after calculation of a0
double fpoly(double v, double *a, int N)
{
    double Ps[N]; double sum=0;
    
    //if(pL->CorL)
    //gsl_sf_legendre_Pl_array(N, 2*v/VMAX - 1, Ps);
    //else
    chebyshev_array(N, 2*v/VMAX - 1, Ps);
    
    for( int i=0; i <= N; i++)
        sum -= a[i]*Ps[i];

    return v*v*exp( sum );
}

class velocityIntegralNorm
{
public:
    double *a;
    int N;
    double operator()(double v) const
    {
        return fpoly( v, a, N);
    }
};

class velocityIntegralonV
{
public:
    double *a;
        int N;
    double operator()(double v) const
    {
        return fpoly( v, a, N)/v;
    }
};

class velocityIntegralvSq
{
public:
    double *a;    
    int N;
    double operator()(double v) const
    {
        return fpoly( v, a, N)*v;
    }
};

void initG(int index, double *a)
{

        //velocityIntegralonV vInt;
        //vInt.a = a;
        //vInt.N = index;
        
        double v[1000];
        double FonVvals[1000];
        double FxVvals[1000];
        for(int i=0; i<1000; i++)
        {
            v[i] = (double) i/2.9e5;
            FonVvals[i] = fpoly( v[i], a, index)/v[i];
            FxVvals[i] = FonVvals[i]*v[i]*v[i];//DEIntegrator<velocityIntegralonV>::Integrate(vInt,vmin[i],VMAX,1e-6);
            //std::cout << v[i] << " " << FxVvals[i] << " " << FonVvals[i] << "\n";
        }
        gsl_spline_init(FonVspline,v,FonVvals,1000);
        gsl_spline_init(FxVspline,v,FxVvals,1000);
            
}

double G(double v0, double ve, double vesc, double vmin, int index, double *a)
{
    
    if ( vmin > VMAX )
        return 0;
    
    if ( index == 0 )
        return ( erf( (ve - vmin) / v0 ) + erf( (ve + vmin) / v0 ) ) / (2*ve);
         
    if ( index == 1 )
    {
        double x = 0;
        
        if(-(ve/v0) + vesc/v0 - vmin/v0 > 0) 
            x += v0*((-4*ve*(-3 + pow(ve,2)/pow(v0,2) - (3*pow(vesc,2))/pow(v0,2) + (3*pow(vmin,2))/pow(v0,2)))/v0 + 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*(-erf(ve/v0 - vmin/v0) - erf(ve/v0 + vmin/v0)))/    (2.*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0))); 
        if(ve/v0 + vesc/v0 - vmin/v0 > 0 && ve/v0 - vesc/v0 + vmin/v0 > 0)
            x += (v0*(2*(ve/v0 + vesc/v0 - vmin/v0)*(3 + (ve/v0 + vesc/v0 - vmin/v0)*(-(ve/v0) + (2*vesc)/v0 + vmin/v0)) + 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*(-erf(vesc/v0) - erf(ve/v0 - vmin/v0))))/(2.*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0)));
    
        return x/v0;
    }
    
           
    if ( index > 1 && index < 8 )
    {
        //velocityIntegralonV vInt;
        //vInt.a = a;
        //vInt.N = index;
        
        if(a[0]==0)
        {
            velocityIntegralNorm vIntN;
            vIntN.a = a;
            vIntN.N = index;
            double norm = DEIntegrator<velocityIntegralNorm>::Integrate(vIntN,0,VMAX,1e-6);
            a[0] = log(norm);
            initG(index,a);
        }
        return gsl_spline_eval_integ(FonVspline, vmin, VMAX, FonVaccel);
        
        //  return DEIntegrator<velocityIntegralonV>::Integrate(vInt,vmin,VMAX,1e-6);
    }
    else        
    {
        std::cout << "vDist not found: " << index << std::endl;
        return(0);
    }
}

double Gsq(double v0, double ve, double vesc, double vmin, int index, double *a)
{
     if ( vmin > VMAX )
        return 0;
    
    if ( index == 0 )
        return (-((vmin/v0 - ve/v0)/exp(pow(vmin/v0 + ve/v0,2))) + (vmin/v0 + ve/v0)/exp(pow(-vmin/v0 + ve/v0,2)))/(2.*sqrt(M_PI)*ve/v0) + ((1 + 2*pow(ve/v0,2))*(-erf(vmin/v0 - ve/v0) + erf(vmin/v0 + ve/v0)))/(2.*ve);
        
    if ( index == 1 )
    {
         double x = 0;
        if(-(ve/v0) + vesc/v0 - vmin/v0 > 0) 
            x += -(v0*(-30*exp(pow(vesc,2)/pow(v0,2) + pow(-(ve/v0) + vmin/v0,2))*(-(ve/v0) + vmin/v0) + 30*exp(pow(vesc,2)/pow(v0,2) + pow(ve/v0 + vmin/v0,2))*(ve/v0 + vmin/v0) +  (4*exp(2*(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2)))*ve*(pow(ve,4)/pow(v0,4) - (10*pow(ve,2)*(1 + pow(vesc,2)/pow(v0,2)))/pow(v0,2) - 15*(2 + (2*pow(vesc,2))/pow(v0,2) + pow(vesc,4)/pow(v0,4) - pow(vmin,4)/pow(v0,4))))/v0 + 15*exp(pow(vesc,2)/pow(v0,2) + 2*(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2)))*sqrt(M_PI)*(1 + (2*pow(ve,2))/pow(v0,2))*(erf(ve/v0 - vmin/v0) + erf(ve/v0 + vmin/v0))))/(20.*exp(2*(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2)))*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0)));
        if(ve/v0 + vesc/v0 - vmin/v0 > 0 && ve/v0 - vesc/v0 + vmin/v0 > 0)
            x +=  (exp(-(pow(ve,2)/pow(v0,2)) - pow(vmin,2)/pow(v0,2))*v0*(2*exp(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2))*(-(pow(ve,5)/pow(v0,5)) + (15*vesc)/v0 + (10*pow(vesc,3))/pow(v0,3) + (4*pow(vesc,5))/pow(v0,5) + (10*pow(ve,3)*(1 + pow(vesc,2)/pow(v0,2)))/pow(v0,3) - (10*pow(vmin,3))/pow(v0,3) - (10*pow(vesc,2)*pow(vmin,3))/pow(v0,5) + (6*pow(vmin,5))/pow(v0,5) + (10*pow(ve,2)*((3*vesc)/v0 + (2*pow(vesc,3))/pow(v0,3) + pow(vmin,3)/pow(v0,3)))/pow(v0,2) + (15*ve*(2 + (2*pow(vesc,2))/pow(v0,2) + pow(vesc,4)/pow(v0,4) - pow(vmin,4)/pow(v0,4)))/v0) + 15*exp(pow(vesc,2)/pow(v0,2))*(-2*exp((2*ve*vmin)/pow(v0,2))*(ve/v0 + vmin/v0) + exp(pow(ve,2)/pow(v0,2) + pow(vmin,2)/pow(v0,2))*sqrt(M_PI)*(1 + (2*pow(ve,2))/pow(v0,2))*(-erf(vesc/v0) - erf(ve/v0 - vmin/v0)))))/(20.*ve*((6*vesc)/v0 + (4*pow(vesc,3))/pow(v0,3) - 3*exp(pow(vesc,2)/pow(v0,2))*sqrt(M_PI)*erf(vesc/v0)));
    
        return x*v0;
     }
        
    if ( index > 1 && index < 8 )
    {
        //velocityIntegralvSq vInt;
        //vInt.a = a;
        //vInt.N = index;
        
        if(a[0]==0)
        {
            velocityIntegralNorm vIntN;
            vIntN.a = a;
            vIntN.N = index;
            double norm = DEIntegrator<velocityIntegralNorm>::Integrate(vIntN,0,VMAX,1e-6);
            a[0] = log(norm);
            initG(index,a);
        }
        return gsl_spline_eval_integ(FxVspline, vmin, VMAX, FxVaccel);
        
        //return DEIntegrator<velocityIntegralvSq>::Integrate(vInt,vmin,VMAX,1e-6);
    }
    else
    {
        std::cout << "vDist not found" << std::endl;
        return(0);
    }

}
