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
        parameterList pL = {};
        char temp[1000];                                   //temporary buffer array
        int myrank = 0;						   			   //MPI process rank
         
        void *pointer;                                     //a pointer for passing parameterList to Multinest
        pointer = (void *) &pL;
     
        //get sampling parameters from file
        int mode = getSamplingPars( &pL, filename); 
        pL.printPars();
        if ( mode < 0 ) 
        {
            std::cout << "Problem with config file, aborting" << std::endl;
            return 0;
        }
        
        // set some MultiNest sampling parameters
        int pWrap[pL.p.nPar];              // which parameters to have periodic boundary conditions?
        int seed = -1;			   	            // random no. generator seed, if < 0 then take the seed from system clock
        int simSeed;
        int ndims = pL.p.nPar;                                // number of parameters in the reconstruction
        int npar  = pL.p.nPar;                             // global number of WIMP parameters
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
	            simSeed = (int) 60000*gsl_rng_uniform(r);		//this is to make sure each MPI thread has the same 'random' data (it would be nice to have more than 60,000 different experiments.. the limit is because it only takes an integer as a seed)
	        }
	        MPI_Bcast(&simSeed,1,MPI_INT,0,MPI_COMM_WORLD);
	        MPI_Barrier(MPI_COMM_WORLD);
        #endif

    //Parameter reconstruction
    if(mode == 0)
    {	 

        //if resuming from previous job
        if( (int)pL.sampling[4] == 1)    	
        { 
            if (pL.w.asimov==1)
            {
                pL.w.asimov=0; pL.sampling[4]=0;
                if(myrank==0)
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
       	    
       	    if( pL.binlessL )
       	    {
       	        if (pL.w.asimov==0)
       	        {
       	            if(myrank==0)
           	            std::cout << "Cannot (as yet) simulate unbinned Asimov data set, switching to Monte-Carlo.." << std::endl;
       	            pL.w.asimov=1; 
       	        }
       	        if(myrank==0)
               	    std::cout << "***experimental feature***" << std::endl;
           	    err = generateUnbinnedData( &(pL.w), &(pL.detectors[i]), 1, simSeed);
           	}
           	else
           	    err = generateBinnedData( &(pL.w), &(pL.detectors[i]), 1, simSeed);
           	    
       	    if(err)
       	        return 1;
       	    
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
            err = writeSamplingOutput(pL);
        
        if(err)
            std::cout << "problem writing output file" << std::endl;
        
        return 0;
    }
    
    //Recoil spectrum
    if(mode == 1)
    {
        int sizeData = 100; //change this if you change arrays below
        double signal[100];
        double background[100]; 
        double Er[100];
        
        for(int j=0; j< pL.ndet; j++)
        {
            std::cout << setiosflags(std::ios::scientific) << setprecision(4);
            if(pL.sampling[3])
            {
                std::cout << "recoil spectrum for detector  " << pL.detectors[j].name << std::endl; 
                std::cout << "Er(keV)       WIMP-rate     Bg-rate      total-rate (/keV/t/year)" << std::endl;
            }
          
            for(int i=1;i<100;i++)
            {
                Er[i] = pL.detectors[j].ErL + (double)i*(pL.detectors[j].ErU-pL.detectors[j].ErL)/100;
                signal[i] = diffWIMPrate(Er[i], &(pL.w), &(pL.detectors[j])); 
                background[i] = diffBgRate(pL.detectors[j],Er[i]);
               
                if(pL.sampling[3])
                    std::cout << Er[i] << "    " << signal[i] << "    " << background[i] << "    " << signal[i]+background[i] << std::endl;
            }
                       
            err = writeRateOutput(pL,j,Er,signal,background,sizeData);
        
            if(err)
                std::cout << "problem writing output file" << std::endl;
        
        }
        return 0;
    }
    
    if( mode!=0 && mode!=1 )
    {
            std::cout << "error with mode parameter, check parameter file" << std::endl;
            return 1;
    }
    
    return 0;
}

