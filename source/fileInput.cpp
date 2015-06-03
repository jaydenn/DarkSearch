#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include <assert.h>
#include "parameterStruct.h"
#include "detectors.h"

using namespace std;

void initialiseArray(char *filename);

//gets sampling parameters from file
int getSamplingPars(parameterList *pL, char *filename) 
{

    FILE* input;
    input = fopen(filename,"r");
    if(input==NULL) 
    {
        printf("unable to open parameter file: %s\n",filename);
        return -1;	
    }
  
    int mode;
    char temp[400];
    char root[50];
    fgets(temp,200,input);

    //Mode switch
    fgets(temp,200,input);
    sscanf(temp,"%d",&mode);
    
    // root for output files
    fgets(temp,200,input);
    sscanf(temp,"%s %*s",root);
    sprintf(pL->root, "%s" , root);		 
    
    //Multinest sampling parameters
    fgets(temp,200,input);
    for (int i=0;i<8;i++)
    {
        fgets(temp,200,input);
        sscanf(temp,"%lf",&(pL->sampling[i]));
    }
  
    //use a binless likelihood
        fgets(temp,200,input);
        sscanf(temp,"%d",&(pL->binlessL));
    //include neutrino background
        fgets(temp,200,input);
        fgets(temp,200,input);
        sscanf(temp,"%d",&(pL->p.nuBg));
    
    //search ranges and nuisance parameters
    char prior[10];
    fgets(temp,200,input);
    fgets(temp,200,input);
    int N = pL->p.nPar; 
    int i = pL->p.nPar; 
    int ind = 1;
    
    //this mess of code makes sure that each parameter can be read in any order and multinest will only loop over the ones with a prior, the others will still be output to file
    fgets(temp,200,input);
    
    while(temp[0]!='/')
    {
       
        if(temp[0]=='M')
        {
            pL->p.Mx[3]= (double)N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.Mx[0]),&(pL->p.Mx[1]),prior);
            if(strcmp(prior,"log")==0) pL->p.Mx[2]=0;
            else if(strcmp(prior,"linear")==0) pL->p.Mx[2]=1;
            else if(strcmp(prior,"gaussian")==0) pL->p.Mx[2]=2;
            else if(strcmp(prior,"none")==0) { pL->p.Mx[2]=3; pL->p.nPar--; i++; pL->p.Mx[3]= pL->p.nPar;}
            else {printf("invalid prior type for Mx\n"); assert(0);}
            sprintf(pL->p.parNames[(int)pL->p.Mx[3]],"Mx");
        }
        
        
        while(temp[0]=='C')
        {        
            pL->p.coeff[ind][3]= (double)N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeff[ind][0]),&(pL->p.coeff[ind][1]),prior);
            if(strcmp(prior,"log")==0) pL->p.coeff[ind][2]=0;
            else if(strcmp(prior,"linear")==0) pL->p.coeff[ind][2]=1;
            else if(strcmp(prior,"gaussian")==0) pL->p.coeff[ind][2]=2; 
            else if(strcmp(prior,"none")==0) { pL->p.coeff[ind][2]=3; pL->p.nPar--; i++; pL->p.coeff[ind][3]= (pL->p.nPar);}
            else {printf("invalid prior type for V%d\n",ind); assert(0);}
            sprintf(pL->p.parNames[(int)pL->p.coeff[ind][3]],"C%d",ind);
            fgets(temp,200,input);
            ind++;
        }
        
     
        if(temp[0]=='r')
        {
            pL->p.rho[3]= N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.rho[0]),&(pL->p.rho[1]),prior);
            pL->p.vSp[0]=pL->p.vSp[0]/3e5;
            pL->p.vSp[1]=pL->p.vSp[1]/3e5;
            if(prior[2]=='g') pL->p.rho[2]=0;
            else if(prior[1]=='i') pL->p.rho[2]=1;
            else if(prior[0]=='g') pL->p.rho[2]=2;
            else if(prior[0]=='n') { pL->p.rho[2]=3; pL->p.nPar--; i++; pL->p.rho[3]= pL->p.nPar;}
            else {printf("invalid prior type for rho\n"); assert(0);}
            sprintf(pL->p.parNames[(int)pL->p.rho[3]],"rho");
        }

        if(temp[1]=='0')
        {
            pL->p.v0[3]= N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.v0[0]),&(pL->p.v0[1]),prior);
            pL->p.v0[0]=pL->p.v0[0]/3e5;
            pL->p.v0[1]=pL->p.v0[1]/3e5;
            if(prior[2]=='g') pL->p.v0[2]=0;
            else if(prior[1]=='i') pL->p.v0[2]=1;
            else if(prior[0]=='g') pL->p.v0[2]=2;
            else if(prior[0]=='n') { pL->p.v0[2]=3; pL->p.nPar--; i++; pL->p.v0[3]= pL->p.nPar;}
            else {printf("invalid prior type for v0\n"); assert(0);}
            sprintf(pL->p.parNames[(int)pL->p.v0[3]],"v0");
        }
        
        if(temp[3]=='c')
        {
            pL->p.vesc[3]= N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.vesc[0]),&(pL->p.vesc[1]),prior);
            pL->p.vesc[0]=pL->p.vesc[0]/3e5;
            pL->p.vesc[1]=pL->p.vesc[1]/3e5;
            if(prior[2]=='g') pL->p.vesc[2]=0;
            else if(prior[1]=='i') pL->p.vesc[2]=1;
            else if(prior[0]=='g') pL->p.vesc[2]=2; 
            else if(prior[0]=='n') { pL->p.vesc[2]=3; pL->p.nPar--; i++; pL->p.vesc[3]= pL->p.nPar;}
            else {printf("invalid prior type for vesc\n"); assert(0);}
            sprintf(pL->p.parNames[(int)pL->p.vesc[3]],"vesc");
        }
        
        if(temp[0]=='v'&&temp[1]=='S')
        {
            pL->p.vSp[3]= N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.vSp[0]),&(pL->p.vSp[1]),prior);
            pL->p.vSp[0]=pL->p.vSp[0]/3e5;
            pL->p.vSp[1]=pL->p.vSp[1]/3e5;
            if(prior[2]=='g') pL->p.vSp[2]=0;
            else if(prior[1]=='i') pL->p.vSp[2]=1;
            else if(prior[0]=='g') pL->p.vSp[2]=2; 
            else if(prior[0]=='n') { pL->p.vSp[2]=3; pL->p.nPar--; i++; pL->p.vSp[3]= pL->p.nPar;}
            else {printf("invalid prior type for vSp\n"); assert(0);}
            sprintf(pL->p.parNames[(int)pL->p.vSp[3]],"vSp");
        }
        
        if(temp[0]=='v'&&temp[1]=='E')
        {
            pL->p.vEp[3]= N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.vEp[0]),&(pL->p.vEp[1]),prior);
            pL->p.vEp[0]=pL->p.vEp[0]/3e5;
            pL->p.vEp[1]=pL->p.vEp[1]/3e5;
            if(prior[2]=='g') pL->p.vEp[2]=0;
            else if(prior[1]=='i') pL->p.vEp[2]=1;
            else if(prior[0]=='g') pL->p.vEp[2]=2; 
            else if(prior[0]=='n') { pL->p.vEp[2]=3; pL->p.nPar--; i++; pL->p.vEp[3]= pL->p.nPar;}
            else {printf("invalid prior type for vEp\n"); assert(0);}
            sprintf(pL->p.parNames[(int)pL->p.vEp[3]],"vEp");
        }
 
        fgets(temp,200,input);
    }
        
    //Detector setup
    char name[20];
    double exp;
    fgets(temp,200,input);
    
    while(temp[0]=='#' && pL->ndet<10)
    {   
        sscanf(temp,"%s %lf", name, &exp);
        newDetector( &(pL->detectors[pL->ndet]), name, exp, pL->ndet);
        pL->ndet++;
        fgets(temp,200,input);
    }    
    
    if (pL->ndet==10) printf("Maximum 10 detectors allowed\n");  
    
    //WIMP parameters
        fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.Mx));
        fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.spin));
        fgets(temp,200,input);
        ind=1;
        while(temp[0]=='C')    
        {    
            sscanf(temp,"%*s %lf",&(pL->w.coeff[ind]));
            fgets(temp,200,input);
            ind++;
        }
        fgets(temp,200,input);
        sscanf(temp,"%d",&(pL->w.asimov));
    
    fclose(input); 

    return mode;
}
