#include <cmath>

//nuclear form factor as a function of recoil energy
double ffactorSI(double A, double Er, int type) 
{

    //Standard Helm 
       
    return 9.05322e6 * exp(-1.9394e-5 * A * Er ) * (-0.00692 * sqrt(2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2)) * sqrt(A * Er) * cos(0.00692 * sqrt(2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2)) * sqrt(A * Er)) + sin( 0.00692 * sqrt(2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2)) * sqrt(A * Er) )) / ( pow( 2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2),1.5) * pow(A * Er,1.5));
        
    //Helm with DarkSUSY implementation
    /*	
        double q = sqrt( 2 * 0.93146 * A * Er / 1e6 );
        double r2 = pow( (1.23*exp(log(A) / 3) - 0.6),2) + 6.22706;
        double x;
        if ( (r2 - 4.05) > 0 )
        {
            x = abs(q * sqrt(r2 - 4.05) / .197327);
        }
        else
        {
            x = abs(q * sqrt(r2) / .197327);
        }
        if ( x > 5e-8)
        {
            return 3 * ( sin(x) - x * cos(x) ) / pow(x,3) * exp( - pow(4.56096 * q,2) / 2 );
        }
        else
        {
            double x2 = pow(x,2);
            return 1 + x2 * ( -1 + x2 * ( 0.00357143 + x2 * ( -0.0000661376 + x2 / 1330560 )));
        }
    */

}
