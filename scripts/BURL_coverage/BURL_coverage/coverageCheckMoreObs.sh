#!/bin/bash

build/burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_fullmodel_moreObs_" \
    -c 10000000 \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 100 \
    -p T \
    -i T

build/burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_withoutphylo_moreObs_" \
    -c 10000000 \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 100 \
    -p F \
    -i T

build/burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_withoutintra_moreObs_" \
    -c 10000000 \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 100 \
    -p T \
    -i F
