# Define the compiler and the compile options
cc       = g++
copt     = -O2 -std=c++17 -march=native -DHYDRA_USE_MKL -DMKL_ILP64 -m64 -I"${MKLROOT}/include"

# Set the directory of lila and hydra here
hydradir = /home/awietek/Research/Software/hydra

# On Linux, you can use the lapack and blas libraries
# lapacklib   = -llapack -lblas
lapacklib   = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
# # on MacOs, you can use the Accelerate framework
# lapacklib   = -framework Accelerate

all:
	$(cc) $(copt) main.cpp -o main -I$(hydradir) -L$(hydradir)/lib  $(lapacklib) -lhydra -fopenmp
