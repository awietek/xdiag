#!/bin/bash

pprof --pdf --ignore kmp benchmark benchmark.prof > benchmark.pdf
