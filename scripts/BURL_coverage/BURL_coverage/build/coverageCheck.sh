#!/bin/bash

./burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_fullmodel_" \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 10 \
    -p T \
    -i T

./burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_withoutphylo_" \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 10 \
    -p F \
    -i T

./burlc \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/coverageCheck/coverage_withoutintra" \
    -nreps 100 \
    -ntips 8 \
    -ntraits 8 \
    -nimp 10 \
    -nobs 10 \
    -p T \
    -i F