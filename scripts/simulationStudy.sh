#!/bin/bash

ntips=(8 16 32)
ntraits=(8)
nimp=(0 8 16)
nobsv=(4 8 16)

for t in "${ntips[@]}"
do
    for tr in "${ntraits[@]}"
    do
        for im in "${nimp[@]}"
        do
            for no in "${nobsv[@]}"
            do
                BURL_coverage/BURL_coverage/build/burlc \
                    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/simulation_study/simulated_${t}_tips_${tr}_nimp_${im}_nobs_${no}_traits_20_reps" \
                    -nreps 20 \
                    -ntips "$t" \
                    -ntraits "$tr" \
                    -nimp "$im" \
                    -nobs "$no" \
                    -p T \
                    -i T
            done
        done
    done
done
