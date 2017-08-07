//definition of parameter structs
#include "detectorStruct.h"
#define PARAMETERSTRUCT_H

struct reconstructionParameters {
    //each parameter array contains:
    // {min, max, prior, parIndex};
    double Mx[4];
    double spin[4];
    double coeffn[16][4];
    double coeffp[16][4];
    double delta[4];
    double rho[4];
    double v0[4];
    double vesc[4];
    double vSp[4];
    double vEp[4];
    char parNames[50][10];
    int isv;
    int nPar;     //number of WIMP parameters tracked
    int nDim;     //dimension of reconstruction
    int vDindex;
    double vLa[8][4];
    
    void printReconstPars()
    {
        
        printf("%d reconstruction parameters:\n",nDim);
        printf("   %.0f < Mx  < %.0f  prior: %.0f, index: %.0f\n",Mx[0],Mx[1],Mx[2],Mx[3]);
        for (int i=1; i<=15; i++)
        {
            printf("   %.0f < c%d_n  < %.0f  prior: %.0f, index: %.0f\n",coeffn[i][0],i,coeffn[i][1],coeffn[i][2],coeffn[i][3]);
            printf("   %.0f < c%d_p  < %.0f  prior: %.0f, index: %.0f\n",coeffp[i][0],i,coeffp[i][1],coeffp[i][2],coeffp[i][3]);
        }
        printf("   %.0f < delta  < %.0f  prior: %.0f, index: %.0f\n",delta[0],delta[1],delta[2],delta[3]);
        printf("   rho = %.2f +\\- %.2f  prior: %.0f, index: %.0f\n",rho[0],rho[1],rho[2],rho[3]);        
        printf("   reconstruct using vDist: %d\n",vDindex);
        printf("   v0  = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",v0[0]*3e5,v0[1]*3e5,v0[2],v0[3]);
        printf("   vesc = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",vesc[0]*3e5,vesc[1]*3e5,vesc[2],vesc[3]);
        printf("   vSp  = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",vSp[0]*3e5,vSp[1]*3e5,vSp[2],vSp[3]);
        printf("   vEp = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",vEp[0]*3e5,vEp[1]*3e5,vEp[2],vEp[3]);
        for (int i=0; i<=vDindex; i++)
            printf("   a%d = %.0f +\\- %.0f  prior: %.0f, index: %.0f\n",i,vLa[i][0],vLa[i][1],vLa[i][2],vLa[i][3]);
    }
    
};

struct WIMPpars {
    double Mx;
    double spin; 
    double coeffp[16];
    double coeffn[16];
    double delta;
    double rho;
    double v0;
    double vesc;
    double vSp;
    double vEp;
    int vDindex;
    double vLa[8];
    int asimov;
    
    void printWIMPpars()
    {
        printf("WIMP parameters:\n");
        printf("   Mx  = %.0f\n",Mx);
        printf("   spin= %.1f\n",spin);
        printf("   delta=%.1f\n",delta);
        printf("   rho = %.2f\n",rho);        
        printf("   v0  = %E\n",v0);
        printf("   vesc= %E\n",vesc);
        printf("   vSp= %E\n",vSp);
        printf("   vEp= %E\n",vEp);
        printf("   asmv= %d\n",asimov);
        printf("   vDist= %d\n",vDindex);
        printf("   a_i={%.0f %.0f %.0f %.0f %.0f}\n",vLa[0],vLa[1],vLa[2],vLa[3],vLa[4]);
        for(int i=1;i<16;i++)
        {
            printf("    c%d=%E   %E\n",i,coeffn[i],coeffp[i]);
        }
    }
    
};

struct parameterList {
    double sampling[11];
    char root[100];
    int binlessL;
    int nbins;
    int ndet;
    
    reconstructionParameters p;
    WIMPpars w;
    detector detectors[10];
    
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
        p.printReconstPars();
        
        printf("Detectors:");
        for(int i=0;i<ndet;i++)
        {    
            printf("%d.",i);
            detectors[i].printDetSpecs();
        }    
    }
        
};
    
