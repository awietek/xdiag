#!/bin/bash

profiler=/home/awietek/Research/Software/gperftools/build/libprofiler.so
nsites=36

LD_PRELOAD=$profiler CPUPROFILE=main.prof ./build/main $nsites

pprof --pdf build/main main.prof > main.prof.pdf
