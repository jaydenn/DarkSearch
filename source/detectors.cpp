#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include "DEIntegrator.h"
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "nuBackground.h"

//returns total # background events per kg/day for bg type, with recoil  Er_min < Er < Er_max
double BgRate(detector det, double Er_min, double Er_max)  						  
{
    if(Er_min < det.ErL)
    { printf("Warning, E_min=%lf < ErL\n",Er_min);}
    
    return gsl_spline_eval_integ(det.background, Er_min, Er_max, det.accel);
}

//only need to calculate background once, store in a table for interpolation, stored as N/kg/day/keV
int InitializeBackground(detector *det)
{
    //get values of bg at relevant energies
    double Er[220];
    double bg[220];
    for(int i=0; i<220; i++)
    { 
        Er[i] = det->ErL + (i-2)*(det->ErU-det->ErL)/200;
        bg[i] = dNnudE(*det, Er[i]) + detBackground(Er[i], det->bg);
    }
    
    //create gsl interpolation object
    return gsl_spline_init(det->background,Er,bg,220);
    
}

int newDetector(detector *det, char *name, double exp, int ndet)
{   
    det->exposure = exp;
    
    FILE *detsINI;
    detsINI = fopen("detectors.ini","r");
    if(detsINI==NULL)
    {
        printf("unable to open detectors.ini\n");
        return -1;	
    }
    
    char temp[200];
    fgets(temp,200,detsINI);
    
    while(strcmp(temp,name)!=0)
    {
        fscanf(detsINI,"%s",temp);
        
        if(feof(detsINI))
        {
            printf("detector '%s' not found\n",name); 
            fclose(detsINI);
            return -1;
        }
    }
    
    sprintf( det->name, "%s", &(name[1]));
    int nIso = det->nIso = 0;
    
    while(temp[0]!='-')
    {
        fscanf(detsINI,"%s",temp);
                     
        if(strcmp(temp,"AM")==0)  
            fscanf(detsINI,"%lf",&(det->AM)); 
        if(strcmp(temp,"Er")==0)  
            fscanf(detsINI,"%lf-%lf",&(det->ErL),&(det->ErU)); 
        if(strcmp(temp,"bg")==0)  
            fscanf(detsINI,"%d",&(det->bg));
        if(strcmp(temp,"eff")==0) 
            fscanf(detsINI,"%d",&(det->eff));         
        if(strcmp(temp,"res")==0) 
            fscanf(detsINI,"%d",&(det->res));
    }
    
    fgets(temp,200,detsINI);
    fgets(temp,200,detsINI);
    fgets(temp,200,detsINI);
    fgets(temp,200,detsINI);
    while(!feof(detsINI) && temp[0]!='_')
    {       
        sscanf(temp,"%d %d %lf",&(det->isoZ[nIso]),&(det->isoA[nIso]),&(det->isoFrac[nIso])); 
        nIso++;
        fgets(temp,200,detsINI);         
    }
    det->nIso = nIso;
             
    fclose(detsINI);
    
    if(InitializeBackground(det))
    { printf("problem initializing background of detector %d\n",ndet); assert(0);}
    
    return 0;
    
}

