#!/bin/bash

pprof --pdf build/main main.prof > benchmark.pdf
