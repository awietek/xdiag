# Define the compiler and the compile options
cc       = clang++
copt     = -O2 -std=c++17

# Set the directory of lila and hydra here
hydradir = /Users/awietek/Research/Software/hydra

# # On Linux, you can use the lapack and blas libraries
# lapacklib   = -llapack -lblas

# On MacOs, you can use the Accelerate framework
lapacklib   = -framework Accelerate

all:
	$(cc) $(copt) -I$(hydradir) -L$(hydradir)/lib  $(lapacklib) -lhydra convert_lfile.cpp -o convert_lfile
