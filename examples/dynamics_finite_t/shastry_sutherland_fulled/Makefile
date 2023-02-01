hydradir = /home/awietek/Research/Software/hydra

ccopt      = -std=c++17 -O2 -march=native -DHYDRA_USE_MKL -DMKL_ILP64 -m64 -I"${MKLROOT}/include" -DHYDRA_USE_HDF5 -lhdf5

lapack =  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
openmp = -fopenmp

all:
	g++ main.cpp -o main $(ccopt) -I$(hydradir) -L$(hydradir)/lib -lhydra $(lapack) $(openmp)


