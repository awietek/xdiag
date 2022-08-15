#!/bin/bash

LD_PRELOAD=/home/awietek/Research/Software/gperftools/build/libprofiler.so CPUPROFILE=benchmark.prof ./benchmark
