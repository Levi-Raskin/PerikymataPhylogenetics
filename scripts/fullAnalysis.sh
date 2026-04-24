#!/bin/bash

### posterior estimation ###
# Lower canine with hominins
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10.tsv" \
    -n 100000000 \
    -p 1000 \
    -s 1000 \
    -c 10
    -nt 10
    -log F

# Lower canine without hominins
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10_no_hominin.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10_no_hominin.tsv" \
    -n 100000000 \
    -p 1000 \
    -s 1000 \
    -c 10
    -nt 10
    -log F

# Lower canine species means
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10_species_means.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10_species_means.tsv" \
    -n 100000000 \
    -p 1000 \
    -s 1000 \
    -c 10
    -nt 10
    -log F


# Upper second incisor
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/ui2/ui2_dec3_10_no_pongo.tsv" \
    -n 100000000 \
    -p 1000 \
    -s 1000 \
    -c 10
    -nt 10
    -log F