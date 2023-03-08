#!/bin/bash
LD_PRELOAD=~/Research/Software/gperftools/build/libprofiler.pc CPUPROFILE=main.prof ./build/main
