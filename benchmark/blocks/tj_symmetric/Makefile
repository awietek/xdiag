arch = flatiron_linux

include ../../../options.mk


all:
	$(cc) $(ccopt) $(lilabackend) -std=c++17 benchmark.cpp -o benchmark $(includes) -I$(liladir) -I../../../ -L../../../lib -lhydra $(libraries)
