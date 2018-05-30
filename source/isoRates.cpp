#include "parameterStruct.h"
#include <iostream>
#include "DMffRates.h"
#include "earthVel.h"

double isoRateT(double Er, WIMPpars *W, int isoA, int isoZ, double T)
{

    double vE = earthVel(T,W->vEp);
    switch(isoZ)
    {
        //FLOURINE
        case 9:
        {
            switch(isoA)
            {
                case 19:
                    return rateZ9A19( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                default:
                    cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
                    return 0;
            }
        }    
        //SODIUM
        case 11:
        {
            switch(isoA)
            {
                case 23:
                    return rateZ11A23( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                default:
                    cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
                    return 0;
            }
        }
        //SILICON
        case 14:
        {
            switch(isoA)
            {
                case 28:
                    return rateZ14A28( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 29:
                    return rateZ14A29( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 30:
                    return rateZ14A30( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                default:
                    cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
                    return 0;
            }
        }
        //ARGON
        //case 18:
        //{   
        //    switch(isoA)
        //    {
        //       case 40:
        //            return rateZ18A40( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
        //        default:
        //            cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
        //            return 0;                
        //   }
        //}    
        //GERMANIUM
        case 32:
        {
            switch(isoA)
            {
                case 70:
                    return rateZ32A70( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 72:
                    return rateZ32A72( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 73:
                    return rateZ32A73( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 74:
                    return rateZ32A74( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                default:
                    cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
                    return 0;
            }
        }
        //IODINE
        case 53:
        {
            switch(isoA)
            {
                case 127:
                    return rateZ53A127( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                default:
                    cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
                    return 0;
            }
        }
        //XENON
        case 54:
        {
            switch(isoA)
            {

                case 129:
                    return rateZ54A129( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 130:
                    return rateZ54A130( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 131:
                    return rateZ54A131( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 132:
                    return rateZ54A132( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 134:
                    return rateZ54A134( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                case 136:
                    return rateZ54A136( Er, W->Mx, W->spin, W->coeffn, W->coeffp, W->rho, W->v0, W->v0+W->vSp+vE, W->vesc, W->vDindex, W->vLa, W->delta);
                default:
                    cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
                    return 0;
            }
        }    
        default:
            cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
            return 0;
    }
        
}

double isoRate(double Er, WIMPpars *W, int isoA, int isoZ)
{
    return isoRateT( Er, W, isoA, isoZ, 0);
}
