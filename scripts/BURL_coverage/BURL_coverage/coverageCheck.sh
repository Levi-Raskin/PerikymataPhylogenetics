#!/bin/bash

build/burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_fullmodel_" \
    -c 10000000 \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 10 \
    -p T \
    -i T &

build/burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_withoutphylo_" \
    -c 10000000 \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 10 \
    -p F \
    -i T &

build/burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_withoutintra_" \
    -c 10000000 \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 10 \
    -p T \
    -i F &

wait
echo "All jobs completed."
