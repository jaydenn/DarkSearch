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
#include "velDist.h"

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
  
    //number of bins to use binless likelihood
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->nbins));
    if (pL->nbins == 0)
        pL->binlessL = 1;
    else
        pL->binlessL = 0;
        
    //include neutrino background
    int nuBg;
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(nuBg));
    
    //isospin violation
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->p.isv));

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
            if(pL->p.isv > 0) //isospin violation allowed
            {
                //neutron coupling will stand in for Cni/Cpi (this allows us to scan log scales while retaining negative interference)
               
                //proton
                sscanf(temp,"%*s %lf %lf %s",&(pL->p.coeffp[ind][0]),&(pL->p.coeffp[ind][1]),prior);
                pL->p.coeffp[ind][2]=priorIndex(prior);
                if(pL->p.coeffp[ind][2]!=3)
                {
                    pL->p.coeffp[ind][3]= (double)pL->p.nDim++;
                    sprintf(pL->p.parNames[(int)pL->p.coeffp[ind][3]],"C%dp",ind);
                    
                    pL->p.coeffn[ind][0]=(double) -pL->p.isv;
                    pL->p.coeffn[ind][1]=(double) +pL->p.isv;
                    pL->p.coeffn[ind][2]=1; //only support linear for now
                    pL->p.coeffn[ind][3]= (double)pL->p.nDim++;
                    sprintf(pL->p.parNames[(int)pL->p.coeffn[ind][3]],"C%dn/C%dp",ind,ind);
                }
                else
                {
                    pL->p.coeffn[ind][2]=3;
                    pL->p.coeffn[ind][3]=-1;
                    pL->p.coeffp[ind][2]=3;
                    pL->p.coeffp[ind][3]=-1;
                }
            }
            else if(pL->p.isv < 0) 
            {
                std::cout << "Incorrect use of isospin parameter, must be zero or positive\n";
                return(-1);
            }
            ret = fgets(temp,200,input);
            ind++;
        }
        
        
        if(temp[0]=='d')
        {
            sscanf(temp,"%*s %lf %lf %s",&(pL->p.delta[0]),&(pL->p.delta[1]),prior);
            pL->p.delta[2]=priorIndex(prior); 
            if(pL->p.delta[2]!=3)
            {
                pL->p.delta[3]= (double)pL->p.nDim++;
                sprintf(pL->p.parNames[(int)pL->p.delta[3]],"delta");
            }
            else
                pL->p.delta[3]=-1;
        
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
            else
                pL->p.rho[3]=-1;
        
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
            else
                pL->p.v0[3]=-1;
                
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
            else
                pL->p.vesc[3]=-1;
            
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
            else
                pL->p.vSp[3]=-1;
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
            else
                pL->p.vEp[3]=-1;
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

                //include non-scanned parameter for a0 (which is fixed by normalization)
                pL->p.vLa[0][3] = (double)pL->p.nDim+pL->p.nPar++;
                pL->p.vLa[0][2] = 3;
                sprintf(pL->p.parNames[(int)pL->p.vLa[0][3]],"a0");
            }
            else
            {
                pL->p.vLa[0][2] = 3;
                pL->p.vLa[0][3] = -1;
            }
        }
        ret = fgets(temp,200,input);
       
    }
    pL->p.nPar+=pL->p.nDim;
    
    //Detector setup
    char name[20];
    double exp;
    ret = fgets(temp,200,input);
    
    while(temp[0]=='#' && pL->ndet<10)
    {
        sscanf(temp,"%s %lf", name, &exp);
        pL->detectors[pL->ndet].nuBg = nuBg;
        if(exp > 0)
        {
            if(newDetector( &(pL->detectors[pL->ndet]), name, exp, pL->ndet))
                return -1;
        }
        else
            std::cout << "detector " << name << " is being ignored\n";
            
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
        sscanf(temp,"%*s %lf",&(pL->w.delta)); 
        ret = fgets(temp,200,input);
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
    /*if(outfile==NULL)
    {
        std::cout << "output file cannot be created" << std::endl;
        return 1;
    }*/
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
    for(int i=0; i<pL.p.nPar; i++)
        outfile  <<  pL.p.parNames[i] << "                         ";

    outfile << std::endl;

    //copy over the Multinest output file
    sprintf(filename,"%s.txt",pL.root);
    infile.open(filename,ios::in);
   /* if(infile==NULL)
    {
        std::cout << "couldn't open Multinest output file" << std::endl;
        return 0;
    }*/
    
    std::string line;
    double prob;
    double totProb = 0;
    double a[8];
    double Gvmin[50][2];
    double gvmin;
    bool first = 1;
    
    while(std::getline(infile,line))
        outfile << line << std::endl;

    
    outfile.close();
    infile.close();
    unlink(filename);

    return 0;
}

int writeVelData(parameterList pL)
{
    char filename[100];
    std::ifstream infile;        
    
    sprintf(filename,"%s.txt",pL.root);
    infile.open(filename,ios::in);
    /*if(infile==NULL)
    {
        std::cout << "couldn't open Multinest output file" << std::endl;
        return 0;
    }*/
    
    std::string line;
    double prob;
    double totProb = 0;
    double a[8];
    double Gvmin[50][2];
    double Fv[50][2];
    double gvmin,fv;
    bool first = 1;
    
    while(std::getline(infile,line))
    {
        sscanf(line.c_str(),"%lf %*lf %*lf %*lf %lf %lf %lf %lf %lf",&prob,&a[1],&a[2],&a[3],&a[4],&a[0]);
        totProb+=prob;
        if( totProb > 0.32 )
        {
            std::cout << totProb << std::endl;
            for (int j=0; j<50; j++)
            {
                gvmin = G(0,0,0,(double)20*j/3e5,4,a);
                fv = fpoly( (double)20*j/3e5, a, 4);
                if (fv < Fv[j][0] || first) 
                    Fv[j][0] = fv;
                if (fv > Fv[j][1] || first) 
                    Fv[j][1] = fv;
                if (gvmin < Gvmin[j][0] || first)
                    Gvmin[j][0] = gvmin;
                if (gvmin > Gvmin[j][1] || first)
                {
                    Gvmin[j][1] = gvmin;
                    first = 0;
                }    
            }
        }
    }
    infile.close();
    
    std::ofstream vOut;
        
    sprintf(filename,"%svelDist.dat",pL.root);
    vOut.open(filename,ios::out);
    for (int j=0; j<50; j++)
        vOut << (double)20*j/3e5 << " " <<  Fv[j][0] << " " <<  Fv[j][1] << " " <<  Gvmin[j][0] << " "  << Gvmin[j][1] << std::endl;
    
    vOut.close();
    
    return 0;
}

/*int writeMonteCarloDataOut(detector *det)
{
    char filename[100];
    std::ofstream outfile;
    
    sprintf(filename,"%s%s_MCevents.dat",pL.root,pL.detectors[detj].name);
    outfile.open(filename,ios::out);
    
    
}*/

int writeRateOutput(parameterList pL, int detj, double *Er, double *signal, double *background, int sizeData)
{
    char filename[100];
    std::ofstream outfile;
    
    sprintf(filename,"%s%s_dRdE.dat",pL.root,pL.detectors[detj].name);
    outfile.open(filename,ios::out);
   if( !outfile )
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
