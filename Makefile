CC = mpif90
<<<<<<< HEAD
FLAGS = -g -Ofast -DMPI -DOMPI_SKIP_MPICXX
=======
FLAGS = -g -O3 -DMPI -DOMPI_SKIP_MPICXX
>>>>>>> 5f981779fc489550def89d0ab81bacb3139ca20b
NESTLIBDIR = ../MultiNest_v3.9

LIBS = -L$(NESTLIBDIR) -lnest3 -llapack -lgsl -lgslcblas -lstdc++ -lgfortran -lm
INCLUDE = -I./source/include 
OBJECTS = source/DarkSearch.o source/nuBackground.o source/detectorFunctions.o source/detectors.o source/directDet.o source/fileIO.o source/formfactorSI.o source/monteCarlo.o source/likelihood.o source/isoRates.o

default: DarkSearch

DarkSearch: $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

source/%.o: source/%.cpp
	$(CC) $< $(FLAGS) $(INCLUDE) -c -o $@

clean:
	-rm source/*.o
	-rm -f ./DarkSearch
