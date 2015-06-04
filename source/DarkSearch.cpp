/*
*      Dark Search v3: Direct Detection Likelihood Calculator
*     -----------------------------------------------------------
*			  Jayden Newstead ASU, Apr. 2015
*
*       A program for the calculation of:
*        - Dark matter direct detection rates and experimental reach
*        - Monte-carlo/Asimov simulation of DM events with model reconstruction
*          with any no. of detectors with generalized nuclear form factors.
*
*/

#include "DarkSearch.h"

int main(int argc, char *argv[])
{	

    char filename[100];
    switch(argc)
    {
        case 1:
        {
            sprintf(filename,"config.dat"); 
            break;
        }
        case 2:
        {
            sprintf(filename,"%s",argv[1]);
            break;
        }
        default:
        {
            std::cout << "Usage: ./DarkSearch {optional: parameter file name}";
            return 0;
        }
    }		  
	 
        //Setup arrays and parameters
        parameterList pL;
        char temp[1000];                                   //temporary buffer array
        int myrank = 0;						   			   //MPI process rank
         
        void *pointer;                                     //a pointer for passing parameterList to Multinest
        pointer = (void *) &pL;
     
        //get sampling parameters from file
        int mode = getSamplingPars( &pL, filename); 
        if ( mode < 0 ) 
        {
            std::cout << "Problem with config file, aborting" << std::endl;
            return 0;
        }
        
        // set some MultiNest sampling parameters
        int pWrap[] = {0,0,0,0,0,0,0,0,0,0,0,0};              // which parameters to have periodic boundary conditions?
        int seed = -1;			   					          // random no. generator seed, if < 0 then take the seed from system clock
        
        int ndims = pL.p.nPar;                                // any combination of mass, sigmaSI, sigmaSIvec, sigmaSD, delta, fn/fp, bn/bp, an/bp and rho_DM, v0, vesc 
        int npar = 21;                                        // npar can be greater than ndim if you want to get other values from the loglike function output to file
        double logZero = -DBL_MAX;							  // points with loglike < logZero will be ignored by MultiNest
        int initMPI = 0;								      // initialize MPI routines?, relevant only if compiling with MPI
        int outfile = 1;								      // write output files?
        int updateInt = 1000000;							  // update interval (for calling dumper, which isn't used here)
        
        int err = 0;                                          //error flag

        //Initialise MPI
        #ifdef MPI
            MPI_Init(&argc, &argv);
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	    if(myrank==0)
	    {
	        gsl_rng_env_setup();
    		T=gsl_rng_default;
	        r=gsl_rng_alloc(T);
	        gsl_rng_set(r, (int)time(NULL));
	        seed = (int) 60000*gsl_rng_uniform(r);		//would be nice to have more than 60,000 different experiments..
	    }
	    MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);
	    MPI_Barrier(MPI_COMM_WORLD);
        #endif

    //Parameter reconstruction
    if(mode == 0)
    {	
        
        //use maximum likelihood values for astro params in simulation of events
        pL.w.rho = pL.p.rho[0];
        pL.w.v0  = pL.p.v0[0];
        pL.w.vesc= pL.p.vesc[0];
        pL.w.vSp = pL.p.vSp[0];
        pL.w.vEp = pL.p.vEp[0];      
       
        //if resuming from previous job
        if( (int)pL.sampling[4] == 1)    	
        { 
            if (pL.w.asimov==1)
            {
                pL.w.asimov=0; pL.sampling[4]=0;
                std::cout << "cannot resume random simulation, switching to Asimov" << std::endl;                   
            }
            else
                std::cout << "resuming from previous job\n" << std::endl; 
        } 
       
        //Generate sim data for each detector
	if(myrank==0)
            std::cout << "Using " << pL.ndet << " detector(s):" << std::endl;
	for(int i=0; i<pL.ndet; i++)
        {                
       	    generateBinnedData( pL.w, pL.p, &(pL.detectors[i]), 1);
	    if(myrank==0)
	            std::cout << "  " << pL.detectors[i].name << " (" << pL.detectors[i].exposure << " t.y): " << pL.detectors[i].nEvents << " events" << std::endl;
      	}   

       
        //run multinest sampling
        if(myrank==0) std::cout << "Starting MultiNest sampling..." << std::endl;
        //  nestRun(               mmodal,                ceff,              nlive,            tol,            efr,ndims, nPar,nCdims,  maxModes,    updInt,          nullZ,     root, seed, pWrap,             feedback,                resume,       outfile,        initMPI, logZero,   loglike, dumper, context)
        nested::run((bool)pL.sampling[0],(bool)pL.sampling[1],(int)pL.sampling[2], pL.sampling[6], pL.sampling[5],ndims, npar, ndims,       100, updateInt, pL.sampling[7],  pL.root, seed, pWrap, (bool)pL.sampling[3], (bool)pL.sampling[4], (bool)outfile, (bool)initMPI, logZero, LogLikedN, dumper, pointer);
        
        #ifdef MPI
            MPI_Finalize();
        #endif
        
        //write parameters to file
        if(myrank==0)
            err = writeOutFile(mode, pL);
        
        if(err)
            cout << "problem writing output file" << endl;
        
        return 0;
    }
    
    //Recoil spectrum
    if(mode == 1)
    {
        double rate = 0; 
        //double rateS= 0;
        char filename[30];
        FILE *output;
        
        pL.w.rho = pL.p.rho[0];
        pL.w.v0 = pL.p.v0[0];
        pL.w.vesc = pL.p.vesc[0];
        pL.w.vSp = pL.p.vSp[0];
        pL.w.vEp = pL.p.vEp[0];        
        
        for(int j=0; j< pL.ndet; j++)
        {
 
            sprintf(filename,"%s%s_dRdE.dat",pL.root,pL.detectors[j].name);
            output = fopen(filename,"w");
            
            fprintf(output,"//recoil spectrum for detector %s\n",pL.detectors[j].name);
            //fprintf(output,"//WIMP pars: Mx = %2.1f SIs = %e SIv = %e SD = %e fnfp = %e bnbp = %e anap = %e del = %e\n// Velocity dist: %d\n", pL.w.Mx, pL.w.SIs, pL.w.SIv, pL.w.SD, pL.w.fnfp, pL.w.bnbp, pL.w.anap, pL.w.del, pL.p.velDist);
            fprintf(output,"//Er(keV)   \trate(/keV/kg/day)\n");
            printf("recoil spectrum for detector %s\nEr(keV)  \trate(/keV/t/year)\n",pL.detectors[j].name);
          
//for(int i=(int)pL.d.detSpecs[j]->ErL; i< (int)pL.d.detSpecs[j]->ErU; i++)
            for(int i=1;i<101;i++)
            {
                if(pL.detectors[j].res > 100)                  //To smear, or not to smear?
                {
                    rate = smearedDiffRate((double)i, pL.w, pL.detectors[j], pL.p);
                }
                else
                {
                    //rate = dNnudE(*(pL.d.detSpecs[0]),(double)i/10); 
                    rate = diffRate((double)i, pL.w, pL.detectors[j], pL.p); 
                }
                
                fprintf(output, "%lf \t%E\n", (double)i, rate);
                if(pL.sampling[3])
                    printf("%E  \t\t%E\n", (double)i, rate);
            }
            fclose(output);
        }
        return 1;
    }
    
    if( mode!=0 && mode!=1 )
    {
            std::cout << "error with mode parameter, check parameter file" << std::endl;
            return 0;
    }
    
    return 1;
}
