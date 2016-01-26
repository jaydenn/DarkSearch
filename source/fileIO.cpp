#include <fstream>
#include <iostream>
#include <iomanip>
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
    char *ret;
    char temp[400];
    char root[100];
    ret = fgets(temp,200,input);

    //Mode switch
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&mode);
    
    // root for output files
    ret = fgets(temp,400,input);
    sscanf(temp,"%s %*s",root);
    sprintf(pL->root, "%s" , root);		 
    
    //Multinest sampling parameters
    ret = fgets(temp,200,input);
    for (int i=0;i<8;i++)
    {
        ret = fgets(temp,200,input);
        sscanf(temp,"%lf",&(pL->sampling[i]));
    }
  
    //use a binless likelihood
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->binlessL));
    
    //include neutrino background
    int nuBg;
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(nuBg));
    
    //isospin violation
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->p.isv));
    if(temp[0]=='p')
        pL->p.isv = 2;
    else if(temp[0]=='n')
        pL->p.isv = 3;
    
    //search ranges and reconstruction parameters    
    //this mess of code makes sure that each parameter can be read in any order and multinest will only loop over the ones with a prior, the others will still be output to file
    ret = fgets(temp,200,input);

    char prior[10];
    pL->p.nPar = 36;
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    int N = pL->p.nPar; 
    int i = pL->p.nPar; 
    int ind = 1;
    
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
        
        if(temp[0]=='s')
        {
            sscanf(temp,"%*s %lf %*s",&(pL->p.spin));
        }
        
        while(temp[0]=='C')
        {
            if(pL->p.isv == 1)
            {
                //neutron
                pL->p.coeffn[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffn[ind][0]),&(pL->p.coeffn[ind][1]),prior);
                if(strcmp(prior,"log")==0) pL->p.coeffn[ind][2]=0;
                else if(strcmp(prior,"linear")==0) pL->p.coeffn[ind][2]=1;
                else if(strcmp(prior,"gaussian")==0) pL->p.coeffn[ind][2]=2; 
                else if(strcmp(prior,"none")==0) { pL->p.coeffn[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffn[ind][3]= (pL->p.nPar);}
                else {printf("invalid prior type for C%d\n",ind); assert(0);}
                sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%dn",ind);
                //proton
                pL->p.coeffp[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffp[ind][0]),&(pL->p.coeffp[ind][1]),prior);
                if(strcmp(prior,"log")==0) pL->p.coeffp[ind][2]=0;
                else if(strcmp(prior,"linear")==0) pL->p.coeffp[ind][2]=1;
                else if(strcmp(prior,"gaussian")==0) pL->p.coeffp[ind][2]=2; 
                else if(strcmp(prior,"none")==0) { pL->p.coeffp[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffp[ind][3]= (pL->p.nPar);}
                else {printf("invalid prior type for C%d\n",ind); assert(0);}
                sprintf(pL->p.parNames[(int)pL->p.coeffp[ind][3]],"C%dp",ind);
            }
            else if(pL->p.isv == 0)
            {
                //neutron
                pL->p.coeffn[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffn[ind][0]),&(pL->p.coeffn[ind][1]),prior);
                if(strcmp(prior,"log")==0) pL->p.coeffn[ind][2]=0;
                else if(strcmp(prior,"linear")==0) pL->p.coeffn[ind][2]=1;
                else if(strcmp(prior,"gaussian")==0) pL->p.coeffn[ind][2]=2; 
                else if(strcmp(prior,"none")==0) { pL->p.coeffn[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffn[ind][3]= (pL->p.nPar);}
                else {printf("invalid prior type for C%d\n",ind); assert(0);}
                sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%dn",ind);
                //proton
                pL->p.coeffp[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffp[ind][0]),&(pL->p.coeffp[ind][1]),prior);
                //not scanning proton values 
                pL->p.coeffp[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffp[ind][3]= (pL->p.nPar);
                sprintf(pL->p.parNames[(int)pL->p.coeffp[ind][3]],"C%dp",ind);
            }
            else if(pL->p.isv == 2)
            {
                //proton
                pL->p.coeffp[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffp[ind][0]),&(pL->p.coeffp[ind][1]),prior);
                if(strcmp(prior,"log")==0) pL->p.coeffp[ind][2]=0;
                else if(strcmp(prior,"linear")==0) pL->p.coeffp[ind][2]=1;
                else if(strcmp(prior,"gaussian")==0) pL->p.coeffp[ind][2]=2; 
                else if(strcmp(prior,"none")==0) { pL->p.coeffp[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffp[ind][3]= (pL->p.nPar);}
                else {printf("invalid prior type for C%d\n",ind); assert(0);}
                sprintf(pL->p.parNames[(int)pL->p.coeffp[ind][3]],"C%dp",ind);
                //not scanning neutron
                pL->p.coeffn[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffn[ind][0]),&(pL->p.coeffn[ind][1]),prior);
                pL->p.coeffn[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffn[ind][3]= (pL->p.nPar);
                sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%dn",ind);
            }
            else if(pL->p.isv == 3)
            {
                //neutron
                pL->p.coeffn[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffn[ind][0]),&(pL->p.coeffn[ind][1]),prior);
                if(strcmp(prior,"log")==0) pL->p.coeffn[ind][2]=0;
                else if(strcmp(prior,"linear")==0) pL->p.coeffn[ind][2]=1;
                else if(strcmp(prior,"gaussian")==0) pL->p.coeffn[ind][2]=2; 
                else if(strcmp(prior,"none")==0) { pL->p.coeffn[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffn[ind][3]= (pL->p.nPar);}
                else {printf("invalid prior type for C%d\n",ind); assert(0);}
                sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%dn",ind);
                //proton
                pL->p.coeffp[ind][3]= (double)N-i--; 
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffp[ind][0]),&(pL->p.coeffp[ind][1]),prior);
                //not scanning proton values 
                pL->p.coeffp[ind][2]=3; pL->p.nPar--; i++; pL->p.coeffp[ind][3]= (pL->p.nPar);
                sprintf(pL->p.parNames[(int)pL->p.coeffp[ind][3]],"C%dp",ind);
            }
            ret = fgets(temp,200,input);
            ind++;
        }
        
     
        if(temp[0]=='r')
        {
            pL->p.rho[3]= N-i--; 
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.rho[0]),&(pL->p.rho[1]),prior);
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
 
        ret = fgets(temp,200,input);
       
    }


    //Detector setup
    char name[20];
    double exp;
    ret = fgets(temp,200,input);
    
    while(temp[0]=='#' && pL->ndet<10)
    {
        sscanf(temp,"%s %lf", name, &exp);
        pL->detectors[pL->ndet].nuBg = nuBg;
        
        if(newDetector( &(pL->detectors[pL->ndet]), name, exp, pL->ndet))
            return -1;
            
        pL->ndet++;
        ret = fgets(temp,200,input);
    }    
    
    if (pL->ndet==10) printf("Maximum 10 detectors allowed\n");  
    
    //WIMP simulation parameters
        ret = fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.Mx));
        ret = fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.spin));
        ret = fgets(temp,200,input);
        ret = fgets(temp,200,input);
        ind=1;
        while(temp[0]=='C')    
        {   
            sscanf(temp,"%*s %lf %lf",&(pL->w.coeffp[ind]),&(pL->w.coeffn[ind]));
            ret = fgets(temp,200,input);
            ind++;
        }
        //astro parameters
        sscanf(temp,"%*s %lf",&(pL->w.rho)); 
        ret = fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.v0)); pL->w.v0/=3e5;
        ret = fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.vesc)); pL->w.vesc/=3e5;
        ret = fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.vSp)); pL->w.vSp/=3e5;
        ret = fgets(temp,200,input);
        sscanf(temp,"%*s %lf",&(pL->w.vEp)); pL->w.vEp/=3e5;
        
        //asimov or random sim?
        ret = fgets(temp,200,input);
        sscanf(temp,"%d",&(pL->w.asimov));
    
    fclose(input); 

    return mode;
}

int writeSamplingOutput(parameterList pL)
{
    char filename[100];
    std::ifstream infile;        
    std::ofstream outfile;
        
    sprintf(filename,"%ssim.dat",pL.root);
    outfile.open(filename,ios::out);
    if(outfile==NULL)
    {
        std::cout << "output file cannot be created" << std::endl;
        return 1;
    }
    if(pL.w.asimov) 
        outfile << "// Monte-carlo simulation" << std::endl;            
    else
        outfile << "// Asimov simulation" << std::endl;

    //print WIMP parameters
    outfile << "// WIMP sim. pars: spin-" << pL.w.spin << ", Mx = " << pL.w.Mx << std::endl << "//           ";
    for( int i=1; i<16; i++)
        outfile << "C" << i << "p = " << pL.w.coeffp[i] << ", ";
    outfile << std::endl << "//           ";
    for( int i=1; i<16; i++)
        outfile << "C" << i << "n = " << pL.w.coeffn[i] << ", ";
    outfile << std::endl;
    outfile << "//Astrophysical parameters: rho = " << pL.p.rho[0] << " GeV/cm^3, v0 = " << pL.p.v0[0]*3E5 << " km/s, vesc = " << pL.p.vesc[0]*3E5 << " km/s" << std::endl; 
        
    //print detector parameters
    outfile << "// Detectors used:" << std::endl;
    for(int i=0; i<pL.ndet; i++)
        outfile << "//      " << pL.detectors[i].name << " - " << pL.detectors[i].exposure << " t.y" << ", nEvents = " << pL.detectors[i].nEvents << std::endl;

    //print parameter headings
    outfile << "//  P                           -2Log(L)                    ";
    for(int i=0; i<36; i++)
        outfile  <<  pL.p.parNames[i] << "                         ";

    outfile << std::endl;

    //copy over the Multinest output file
    sprintf(filename,"%s.txt",pL.root);
    infile.open(filename,ios::in);
    if(infile==NULL)
    {
        std::cout << "couldn't open Multinest output file" << std::endl;
        return 0;
    }
    
    std::string line;
    while(std::getline(infile,line))
        outfile << line << std::endl;
    outfile.close();
    infile.close();
    
    return 0;
}

int writeRateOutput(parameterList pL, int detj, double *Er, double *signal, double *background, int sizeData)
{
    char filename[100];
    std::ofstream outfile;
    
    sprintf(filename,"%s%s_dRdE.dat",pL.root,pL.detectors[detj].name);
    outfile.open(filename,ios::out);
    if(outfile==NULL)
    {
        std::cout << "output file cannot be created" << std::endl;
        return 1;
    }
    
    //write out WIMP parameters
    outfile << "//recoil spectrum for detector "  << pL.detectors[detj].name << std::endl;
    outfile << "//WIMP sim. pars: Mx = " << pL.w.Mx << " GeV, spin = " << pL.w.spin << std::endl << "//           ";
    for( int i=1; i<16; i++)
        outfile << "C" << i << "p = " << pL.w.coeffp[i] << ", ";
    outfile << std::endl << "//           ";
    for( int i=1; i<16; i++)
        outfile << "C" << i << "n = " << pL.w.coeffn[i] << ", ";
    outfile << std::endl;
    outfile << "//Astrophysical parameters: rho = " << pL.p.rho[0] << " GeV/cm^3, v0 = " << pL.p.v0[0]*3E5 << " km/s, vesc = " << pL.p.vesc[0]*3E5 << " km/s" << std::endl; 
    outfile << "//Er(keV)         WIMP-rate         Bg-rate           total-rate (/keV/t/year)" << endl;
    
    outfile << std::setiosflags(std::ios::scientific) << std::setprecision(8);
    //write out rate data
    for (int i=1; i < sizeData; i++)
         outfile << Er[i] <<  "    " << signal[i] << "    " << background[i] << "    " << signal[i]+background[i] << std::endl;

    outfile.close();
    
    return 0;
}
