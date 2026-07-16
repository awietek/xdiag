### MEMORY PROFILE ###
- compile xdiag using

cmake -S . -B build -D CMAKE_CXX_COMPILER=icpx -D CMAKE_CXX_FLAGS="-L/home/awietek/Research/Software/gperftools/build -ltcmalloc -g"

- install xdiag
- compile application using

cmake -S . -B build -D CMAKE_CXX_COMPILER=icpx -D CMAKE_CXX_FLAGS="-L/home/awietek/Research/Software/gperftools/build -ltcmalloc -g"

- run using
HEAPPROFILE=heap ./build/main 28

- visualize using
pprof --pdf ./build/main heap.0001.heap > heap.pdf
