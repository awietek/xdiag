# Define the compiler and the compile options
cc       = clang++
copt     = -O3 -Ofast -mavx -std=c++17 -m64 -Wno-return-type-c-linkage

# Set the directory of lila and hydra here
liladir  = /Users/awietek/Research/Software/lila
hydradir = /Users/awietek/Research/Software/hydra

# Choose the lapack backend for lila
lilabackend = -DLILA_USE_ACCELERATE
lapacklib   = -framework Accelerate

all:
	$(cc) $(copt) -I$(liladir) -I$(hydradir) -L$(hydradir)/lib $(lilabackend) $(lapacklib) -lhydra spinhalf_chain.cpp -o spinhalf_chain
