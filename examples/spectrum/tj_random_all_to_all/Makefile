# Define the compiler and the compile options
cc       = clang++
copt     = -O2 -std=c++17

# Set the directory of lila and hydra here
liladir  = /Users/awietek/Research/Software/lila
hydradir = /Users/awietek/Research/Software/hydra

# Choose the lapack backend for lila
lilabackend = -DLILA_USE_ACCELERATE
lapacklib   = -framework Accelerate

all:
	$(cc) $(copt) -I$(liladir) -I$(hydradir) -L$(hydradir)/lib $(lilabackend) $(lapacklib) -lhydra tj_random_all_to_all.cpp -o tj_random_all_to_all
