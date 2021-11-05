#!/bin/bash

LD_PRELOAD=/mnt/home/awietek/gperftools/lib/libprofiler.so CPUPROFILE=benchmark.prof ./benchmark
