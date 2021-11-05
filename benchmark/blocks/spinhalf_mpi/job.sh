#!/bin/bash
#SBATCH -N4 --exclusive

module purge
module load gcc llvm openmpi4 lib/hdf5 intel/mkl slurm

mpirun ./benchmark
