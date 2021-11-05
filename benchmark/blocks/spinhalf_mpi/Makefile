arch = flatiron_linux

include ../../../options.mk


all:
	$(mpicc) $(ccopt) $(ccarch) benchmark.cpp -o benchmark $(includes) -I../../../ -L../../../lib -lhydra_mpi $(libraries)
