#!/bin/bash

ntips=(8 16 32)
ntraits=(8 16 32)
nimp=(0 8 16 32)
nobsv=(2 4 8 16)

max_jobs=8

nc=10000000
nr=100

run_burlc () {
    local t="$1" tr="$2" im="$3" no="$4"
    local outbase="/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/simulation_study/simulated_tips_${t}_traits_${tr}_nimp_${im}_nobs_${no}_reps_${nr}_"
    local outfile="${outbase}CoverageResults.txt"

    if [ -f "$outfile" ]; then
        echo "Skipping (already exists): $outfile"
        return
    fi

    BURL_coverage/BURL_coverage/build/burlc \
        -o "$outbase" \
        -c "$nc" \
        -nreps "$nr" \
        -ntips "$t" \
        -ntraits "$tr" \
        -nimp "$im" \
        -nobs "$no" \
        -p T \
        -i T
}

run_burlc_wo_intra () {
    local t="$1" tr="$2" im="$3" no="$4"
    local outbase="/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/simulation_study/simulated_wo_intra_${t}_traits_${tr}_nimp_${im}_nobs_${no}_reps_${nr}_"
    local outfile="${outbase}CoverageResults.txt"

    if [ -f "$outfile" ]; then
        echo "Skipping (already exists): $outfile"
        return
    fi

    BURL_coverage/BURL_coverage/build/burlc \
        -o "$outbase" \
        -c "$nc" \
        -nreps "$nr" \
        -ntips "$t" \
        -ntraits "$tr" \
        -nimp "$im" \
        -nobs "$no" \
        -p T \
        -i F
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


for t in "${ntips[@]}"; do
    for tr in "${ntraits[@]}"; do
        for im in "${nimp[@]}"; do
            for no in "${nobsv[@]}"; do

                while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]; do
                    wait -n 2>/dev/null || sleep 5
                done

                run_burlc_wo_intra "$t" "$tr" "$im" "$no" &

            done
        done
    done
done


wait