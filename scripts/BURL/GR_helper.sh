#!/bin/bash

for i in 0 1 2 3; do
    ./burl \
        -i /Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo.csv \
        -it /Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt \
        -o /Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/GR_test_chains_longer/chain${i}.tsv \
        -n 1000000000 \
        -p 1000 \
        -s 1000 \
        -c 10 \
        -nt 10 \
        -log F
done