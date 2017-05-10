#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include "parameterStruct.h"
#include "detectors.h"

using namespace std;

void initialiseArray(char *filename);

int priorIndex(char *prior)
{
    if(strcmp(prior,"log")==0)
        return 0;
    if(strcmp(prior,"linear")==0)
        return 1;
    if(strcmp(prior,"gaussian")==0)
        return 2;
    if(strcmp(prior,"none")==0) 
        return 3;
    else 
        std::cout << "invalid prior type\n"; 
    assert(0);
}

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
    
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->p.vDindex));
    if(pL->p.vDindex > 7)
    {
        std::cout << "maximum polynomial degree supported is 7\n";
        pL->p.vDindex = 7;
    }
    
    //search ranges and reconstruction parameters    
    //this mess of code makes sure that each parameter can be read in any order and multinest will only loop over the ones with a prior
    ret = fgets(temp,200,input);

    char prior[10];
    pL->p.nDim = 0;
    pL->p.nPar = 0;
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    int ind = 1;
    
    while(temp[0]!='/')
    {
        if(temp[0]=='M')
        {
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.Mx[0]),&(pL->p.Mx[1]),prior);
            pL->p.Mx[2]=priorIndex(prior); 
            if(pL->p.Mx[2]!=3)
            {
                pL->p.Mx[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.Mx[3]],"Mx");
            }
        }

        if(temp[0]=='s')
        {
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.spin[0]),&(pL->p.spin[1]),prior);
            pL->p.spin[2]=priorIndex(prior); 
            if(pL->p.spin[2]!=3)
            {
                pL->p.spin[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.spin[3]],"spin");
            }
        }
       
        while(temp[0]=='C')
        {
            if(pL->p.isv == 0) //no isospin violation
            {
                //neutron coupling will stand in for both
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffn[ind][0]),&(pL->p.coeffn[ind][1]),prior);
                pL->p.coeffn[ind][2]=priorIndex(prior);
                if(pL->p.coeffn[ind][2]!=3)
                {
                    pL->p.coeffn[ind][3]= (double)pL->p.nDim++;
                    sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%d",ind);
                }
                //proton set to match neutron
                pL->p.coeffp[ind][0]=pL->p.coeffn[ind][0];
                pL->p.coeffp[ind][1]=pL->p.coeffn[ind][1];
                pL->p.coeffp[ind][2]=3;
            }
            else if(pL->p.isv == 1) //isospin violation allowed
            {
                //neutron
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffn[ind][0]),&(pL->p.coeffn[ind][1]),prior);
                pL->p.coeffn[ind][2]=priorIndex(prior);
                if(pL->p.coeffn[ind][2]!=3)
                {
                    pL->p.coeffn[ind][3]= (double)pL->p.nDim++;
                    sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%dn",ind);
                }
                
                //proton
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffp[ind][0]),&(pL->p.coeffp[ind][1]),prior);
                pL->p.coeffp[ind][2]=priorIndex(prior);
                if(pL->p.coeffp[ind][2]!=3)
                {
                    pL->p.coeffp[ind][3]= (double)pL->p.nDim++;
                    sprintf(pL->p.parNames[(int)pL->p.coeffp[ind][3]],"C%dp",ind);
                }
            }
            else if(pL->p.isv == 2) //proton only
            {
                //proton
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffp[ind][0]),&(pL->p.coeffp[ind][1]),prior);
                pL->p.coeffp[ind][2]=priorIndex(prior);
                if(pL->p.coeffp[ind][2]!=3)
                {
                    pL->p.coeffp[ind][3]= (double)pL->p.nDim++;
                    sprintf(pL->p.parNames[(int)pL->p.coeffp[ind][3]],"C%dp",ind);
                }
                //not scanning neutron
                pL->p.coeffn[ind][0]=0;
                pL->p.coeffn[ind][1]=0;
                pL->p.coeffn[ind][2]=3;
            }
            else if(pL->p.isv == 3) //neutron only
            {
                //neutron
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffn[ind][0]),&(pL->p.coeffn[ind][1]),prior);
                pL->p.coeffn[ind][2]=priorIndex(prior);
                if(pL->p.coeffn[ind][2]!=3)
                {
                    pL->p.coeffn[ind][3]= (double)pL->p.nDim++;
                    sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%dn",ind);
                }
                //proton
                pL->p.coeffp[ind][0]=0;
                pL->p.coeffp[ind][1]=0;
                pL->p.coeffp[ind][2]=3;
            }
            ret = fgets(temp,200,input);
            ind++;
            
        }
        
     
        if(temp[0]=='r')
        {
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.rho[0]),&(pL->p.rho[1]),prior);
            pL->p.rho[2]=priorIndex(prior); 
            if(pL->p.rho[2]!=3)
            {
                pL->p.rho[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.rho[3]],"rho");
            }
        
        }

        if(temp[0]=='v' && temp[1]=='0')
        {
        
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.v0[0]),&(pL->p.v0[1]),prior);
            pL->p.v0[2]=priorIndex(prior); 
            if(pL->p.v0[2]!=3)
            {
                pL->p.v0[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.v0[3]],"v0");
            }
            pL->p.v0[0]=pL->p.v0[0]/3e5;
            pL->p.v0[1]=pL->p.v0[1]/3e5;
            
        }
        
        if(temp[3]=='c')
        {
              
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.vesc[0]),&(pL->p.vesc[1]),prior);
            pL->p.vesc[2]=priorIndex(prior); 
            if(pL->p.vesc[2]!=3)
            {
                pL->p.vesc[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.vesc[3]],"vesc");
            }
            pL->p.vesc[0]=pL->p.vesc[0]/3e5;
            pL->p.vesc[1]=pL->p.vesc[1]/3e5;
            
        }
        
        if(temp[0]=='v'&&temp[1]=='S')
        {
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.vSp[0]),&(pL->p.vSp[1]),prior);
            pL->p.vSp[2]=priorIndex(prior); 
            if(pL->p.vSp[2]!=3)
            {
                pL->p.vSp[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.vSp[3]],"vSp");
            }
            pL->p.vSp[0]=pL->p.vSp[0]/3e5;
            pL->p.vSp[1]=pL->p.vSp[1]/3e5;
        }
        
        if(temp[0]=='v'&&temp[1]=='E')
        {
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.vEp[0]),&(pL->p.vEp[1]),prior);
            pL->p.vEp[2]=priorIndex(prior); 
            if(pL->p.vEp[2]!=3)
            {
                pL->p.vEp[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.vEp[3]],"vSp");
            }
            pL->p.vEp[0]=pL->p.vEp[0]/3e5;
            pL->p.vEp[1]=pL->p.vEp[1]/3e5;
        }
        
        if(temp[0]=='a')
        {
            
            //set same for all a's
            int aInd=1;
            if( pL->p.vDindex > 1 )
            {
                while (aInd <= pL->p.vDindex)
                {
                    sscanf(temp,"%*s %lf %lf %s",&(pL->p.vLa[aInd][0]),&(pL->p.vLa[aInd][1]),prior);
                    pL->p.vLa[aInd][2]=priorIndex(prior); 
                    if(pL->p.vLa[aInd][2]!=3)
                    {
                        pL->p.vLa[aInd][3] = (double)pL->p.nDim++;
                        sprintf(pL->p.parNames[(int)pL->p.vLa[aInd][3]],"a%d",aInd);
                    }
                    aInd++;
                }
            }
            //include non-scanned parameter for a0 (which is fixed by normalization)
            pL->p.vLa[aInd][3] = (double)pL->p.nDim++; pL->p.nPar++;
            pL->p.vLa[aInd][2] = 3;
            sprintf(pL->p.parNames[(int)pL->p.vLa[0][3]],"a0");
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
        ret = fgets(temp,200,input);
        sscanf(temp,"%*s %d",&(pL->w.vDindex));
        
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
    for(int i=0; i<37; i++)
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
    
    unlink(filename);
    
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
