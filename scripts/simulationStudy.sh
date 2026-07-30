#!/bin/bash

ntips=(8 16 32)
ntraits=(8 16 32)
nimp=(0 8 16)
nobsv=(2 4 8 16)

max_jobs=3

run_burlc () {
    local t="$1" tr="$2" im="$3" no="$4"
    BURL_coverage/BURL_coverage/build/burlc \
        -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/simulation_study/simulated_${t}_tips_${tr}_nimp_${im}_nobs_${no}_traits_100_reps" \
        -c 1000000 \
        -nreps 100 \
        -ntips "$t" \
        -ntraits "$tr" \
        -nimp "$im" \
        -nobs "$no" \
        -p T \
        -i T
}

for t in "${ntips[@]}"; do
    for tr in "${ntraits[@]}"; do
        for im in "${nimp[@]}"; do
            for no in "${nobsv[@]}"; do

                while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]; do
                    wait -n 2>/dev/null || sleep 5
                done

                run_burlc "$t" "$tr" "$im" "$no" &

            done
        done
    done
done

wait