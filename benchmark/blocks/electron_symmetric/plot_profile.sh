#!/bin/bash

pprof --pdf benchmark --ignore kmp benchmark.prof > benchmark.pdf
