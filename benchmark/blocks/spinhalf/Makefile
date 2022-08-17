arch = flatiron_linux

include ../../../options.mk

openmp=-DHYDRA_ENABLE_OPENMP -fopenmp

all:
	$(cc) $(ccopt) $(lilabackend) $(openmp) -std=c++17  benchmark.cpp -o benchmark $(includes) -I../../../ -L../../../lib -lhydra $(libraries)
