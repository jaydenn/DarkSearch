#include <math.h>

//returns earth velocity in galactic frame at time T (as a fraction of the year), returns
double earthVel(double T, double vE, double *vS, double v0)
{
    double nn=365.25*T+3652.5; //this needs fixing but is of minor importance, there isn't much T dependence year to year.

    return sqrt(pow(fabs(vS[0] - vE*((0.993821 + 3.603011635865845e-8*nn)*
           (cos(4.894950420143296 + 0.017202792393721557*nn) + 
             0.01671*cos(4.851805881033997 + 0.034404762737365424*nn)) + 
          (0.054857 - 6.634360027378507e-7*nn)*
           (sin(4.894950420143296 + 0.017202792393721557*nn) + 
             0.01671*sin(4.851805881033997 + 0.034404762737365424*nn)))),2) + 
    pow(fabs(v0 + vS[1] - vE*((0.110992 - 3.2446269678302534e-7*nn)*
           (cos(4.894950420143296 + 0.017202792393721557*nn) + 
             0.01671*cos(4.851805881033997 + 0.034404762737365424*nn)) + 
          (-0.494109 - 7.36208076659822e-8*nn)*
           (sin(4.894950420143296 + 0.017202792393721557*nn) + 
             0.01671*sin(4.851805881033997 + 0.034404762737365424*nn)))),2) + 
    pow(fabs(vS[2] - vE*((0.000352 + 5.82258726899384e-7*nn)*
           (cos(4.894950420143296 + 0.017202792393721557*nn) + 
             0.01671*cos(4.851805881033997 + 0.034404762737365424*nn)) + 
          (0.867666 + 4.232717316906228e-11*nn)*
           (sin(4.894950420143296 + 0.017202792393721557*nn) + 
             0.01671*sin(4.851805881033997 + 0.034404762737365424*nn)))),2));
}
