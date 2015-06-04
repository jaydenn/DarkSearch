//definition of parameter structs
#include "detectorStruct.h"
#define PARAMETERSTRUCT_H

struct physicalParameters {
    //each parameter array contains:
    // {min, max, prior, parIndex};
    double Mx[4];
    
    double rho[4];
    double v0[4];
    double vesc[4];
    double vSp[4];
    double vEp[4];
    double coeff[16][4];
    double spin;
    int nPar;
    char parNames[11][5];
    int nuBg;
   
    physicalParameters()
    {
        nPar = 21;
        for(int i=0;i<4;i++)
        {
            Mx[i]=0;
            
            rho[i]=0;
            v0[i]=0;
            vesc[i]=0;
            vSp[i]=0;
            vEp[i]=0;            
        }
        //for(int i=0;i<11;i++)
        //    parNames[i]=T;
        
        nPar=0;
        nuBg=0;
    }
    
    void printPhysPars()
    {
        printf("%d physical parameters:\n",nPar);
        printf("   %.0f < Mx  < %.0f  prior: %.0f, index: %.0f\n",Mx[0],Mx[1],Mx[2],Mx[3]);
        printf("   rho = %.2f +\\- %.2f  prior: %.0f, index: %.0f\n",rho[0],rho[1],rho[2],rho[3]);        
        printf("   v0  = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",v0[0]*3e5,v0[1]*3e5,v0[2],v0[3]);
        printf("   vesc = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",vesc[0]*3e5,vesc[1]*3e5,vesc[2],vesc[3]);
        printf("   vSp  = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",vSp[0]*3e5,vSp[1]*3e5,vSp[2],vSp[3]);
        printf("   vEp = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",vEp[0]*3e5,vEp[1]*3e5,vEp[2],vEp[3]);
        
    }
    
};

struct WIMPpars {
    double Mx;
    double spin; 
    double coeff[16];
    double rho;
    double v0;
    double vesc;
    double vSp;
    double vEp;
    int asimov;
    
    void printWIMPpars()
    {
        printf("WIMP parameters:\n");
        printf("   Mx  = %.0f\n",Mx);
        printf("   spin= %.1f\n",spin);
        printf("   rho = %.2f\n",rho);        
        printf("   v0  = %.0f\n",v0);
        printf("   vesc= %.0f\n",vesc);
        printf("   vSp= %.0f\n",vSp);
        printf("   vEp= %.0f\n",vEp);
        printf("   asmv= %d\n",asimov);
        for(int i=1;i<16;i++)
        {
            printf("    c%d=%lf\n",i,coeff[i]);
        }
    }
    
    WIMPpars()
    {
        rho=0; v0=0; vesc=0; vSp=0; vEp=0;
    }
};

struct parameterList {
    double sampling[11];
    char root[20];
    int binlessL;
    int ndet;
    
    physicalParameters p;
    WIMPpars w;
    detector *detectors;
    
    void printPars()
    {
        printf("root is: %s\n",root);
        printf("sampling:\n");
        printf("   Mode sep   %d\n",(int)sampling[0]);
        printf("   Const eff  %d\n",(int)sampling[1]);
        printf("   live pnts  %d\n",(int)sampling[2]);
        printf("   std out    %d\n",(int)sampling[3]);
        printf("   resume     %d\n",(int)sampling[4]);
        printf("   req. eff   %.1f\n",sampling[5]);
        printf("   tol        %.1f\n",sampling[6]);
        printf("   Ztol       %.2E\n",sampling[7]);
               
        w.printWIMPpars();
        p.printPhysPars();
        
        printf("Detectors:");
        for(int i=0;i<ndet;i++)
        {    printf("%d.",i);
            detectors[i].printDetSpecs();
        }    
    }
    
    parameterList()
    {
        ndet=0;
        p.nPar=21;
        detectors = new detector[10];   
    }
    
};
    
