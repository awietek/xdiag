#!/bin/bash

pprof --pdf --ignore gomp benchmark benchmark.prof > benchmark.pdf
