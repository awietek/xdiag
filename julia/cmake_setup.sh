# Macbook
cmake -S . -B build -D XDIAG_JULIA_WRAPPER=On -D CMAKE_PREFIX_PATH=/Users/awietek/.julia/artifacts/e94dd602bf68f8208260d642c2792a0624e16169 -D CMAKE_CXX_COMPILER=clang-20 -D CMAKE_SHARED_LINKER_FLAGS=-lc++

###!!!!!  export OMP_NUM_THREADS=1
