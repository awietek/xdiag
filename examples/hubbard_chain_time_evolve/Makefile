# Define the compiler and the compile options
cc       = g++
copt     = -O2 -std=c++17

# Set the directory of lila and hydra here
hydradir = /Users/awietek/Research/Software/hydra

# # On Linux, you can use the lapack and blas libraries
# lapacklib   = -llapack -lblas

# on MacOs, you can use the Accelerate framework
lapacklib   = -framework Accelerate

all:
	$(cc) $(copt) main.cpp -o main -I$(hydradir) -L$(hydradir)/lib  $(lapacklib) -lhydra
