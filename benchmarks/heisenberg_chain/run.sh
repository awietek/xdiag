#!/bin/bash

profiler=/home/awietek/Research/Software/gperftools/build/libprofiler.so
nsites=28

LD_PRELOAD=$profiler CPUPROFILE=bench.prof ./build/bench $nsites

pprof --pdf build/bench bench.prof > bench.prof.pdf
