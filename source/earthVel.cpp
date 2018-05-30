#include <math.h>

double earthVel(double T, double vE)
{
    return vE*sin(2*M_PI*T);
}
