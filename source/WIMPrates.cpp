#include "parameterStruct.h"
#include <iostream>
#include "DMffRates.h"

double rate(double Er, WIMPpars W, physicalParameters P, int isoA, int isoZ)
{
    switch(isoZ)
    {
        //FLOURINE
        case 9:
        {
            switch(isoA)
            {
                case 19:
                    return rate9A19( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
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
                    return rate11A23( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
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
                    return rate14A28( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 29:
                    return rate14A29( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 30:
                    return rate14A30( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                default:
                    cout << "isotope A=" << isoA << " Z=" << isoZ << " not found, rate will be inaccurate" << endl;
                    return 0;
            }
        }
        //GERMANIUM
        case 32:
        {
            switch(isoA)
            {
                case 70:
                    return rate32A70( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 72:
                    return rate32A72( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 73:
                    return rate32A73( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 74:
                    return rate32A74( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
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
                    return rate32A74( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
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
                    return rate54A129( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 130:
                    return rate54A130( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 131:
                    return rate54A131( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 132:
                    return rate54A131( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 134:
                    return rate54A131( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
                case 136:
                    return rate54A131( Er, W.Mx, W.spin, W.coeff, W.rho, W.v0, W.vSp, W.vesc);
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

