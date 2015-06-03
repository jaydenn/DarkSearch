#include <math.h>
#include <stdio.h>
double detEff(double Er, int type){
switch( type ) {
case 0: {return 1;}
case 1: {return .4;}
case 2: {return .4;}
case 3: {return .5;}
case 4: {return 1;}
}printf("invalid detetector efficiency\n"); return NAN; }
double detBackground(double Er, int type){
switch( type ) {
case 0: {return 1e-6;}
case 1: {return 1e-9;}
case 2: {return 4e-6;}
case 3: {return 1e-5;}
case 4: {return 1e-6;}
}printf("invalid detetector background\n"); return NAN; }
double detRes(double Er, int type){
switch( type ) {
case 0: {return 5;}
case 1: {return 5;}
case 2: {return 5;}
case 3: {return 5;}
case 4: {return 5;}
}printf("invalid detetector resolution\n"); return NAN; }
