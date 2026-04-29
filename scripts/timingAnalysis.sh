#!/bin/bash

ntips=(8 16 32 64 128 256 512)
ntraits=(10 25 50 100 200)

for t in "${ntips[@]}"
do
    for tr in "${ntraits[@]}"
    do
        BURL_coverage/BURL_coverage/build/burlc \
            -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/timing/timing_${t}_tips_${tr}_traits_10_reps" \
            -c 1000000 \
            -nreps 10 \
            -ntips "$t" \
            -ntraits "$tr" \
            -nimp 10 \
            -nobs 10 \
            -p T \
            -i T
    done
done