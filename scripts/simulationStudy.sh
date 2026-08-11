#!/bin/bash

ntips=(8 16 32)
ntraits=(8 16 32)
nimp=(0 8 16 32)
nobsv=(2 4 8 16)

max_jobs=32

nc=10000000
nr=100

run_burlc () {
    local t="$1" tr="$2" im="$3" no="$4"
    local outbase="/space/s1/levi_raskin/simulation_study/simulated_${t}_tips_${tr}_traits_${im}_nimp_${no}_nobs_${nr}_reps_"
    local outfile="${outbase}CoverageResults.txt"

    if [ -f "$outfile" ]; then
        echo "Skipping (already exists): $outfile"
        return
    fi

    BURL_coverage/build/burlc \
        -o "$outbase" \
        -c "$nc" \
        -nreps "$nr" \
        -ntips "$t" \
        -ntraits "$tr" \
        -nimp "$im" \
        -nobs "$no" \
        -p T \
        -i T \
        -runk F
}

run_burlc_wo_intra () {
    local t="$1" tr="$2" im="$3" no="$4"
    local outbase="/space/s1/levi_raskin/simulation_study/simulated_wo_intra_${t}_tips_${tr}_traits_${im}_nimp_${no}_nobs_${nr}_reps_"
    local outfile="${outbase}CoverageResults.txt"

    if [ -f "$outfile" ]; then
        echo "Skipping (already exists): $outfile"
        return
    fi

    BURL_coverage/build/burlc \
        -o "$outbase" \
        -c "$nc" \
        -nreps "$nr" \
        -ntips "$t" \
        -ntraits "$tr" \
        -nimp "$im" \
        -nobs "$no" \
        -p T \
        -i F \
        -runk F
}

run_burlc_runk () {
    local im="$1" no="$2"
    local outbase="/space/s1/levi_raskin/simulation_study/simulated_runk_8_tips_8_traits_${im}_nimp_${no}_nobs_5000_reps_"
    local outfile="${outbase}CoverageResults.txt"

    if [ -f "$outfile" ]; then
        echo "Skipping (already exists): $outfile"
        return
    fi

    BURL_coverage/build/burlc \
        -o "$outbase" \
        -c "$nc" \
        -nreps 5000 \
        -ntips 8 \
        -ntraits 8 \
        -nimp "$im" \
        -nobs "$no" \
        -p T \
        -i T \
        -runk T
}

### RUNK ###
for im in 0 32; do
    for no in 2 16; do

        while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]; do
            wait -n 2>/dev/null || sleep 5
        done

        run_burlc_runk "$im" "$no" &

    done
done

### Coverage ###
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

### Coverage wo intra ###
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