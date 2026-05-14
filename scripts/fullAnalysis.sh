#!/bin/bash

#number of MCMC cycles
nc=100000000

### posterior estimation ###
# Lower canine with hominins
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10.tsv" \
    -n "$nc" \
    -p 1000 \
    -s 1000 \
    -c 10 \
    -nt 10 \
    -log F

#Gelman-Rubin
for i in 1 2 3 4; do
    BURL/build/burl \
        -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv" \
        -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
        -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/gelmanRubin/out${i}.tsv" \
        -n "$nc" \
        -p 1000 \
        -s 1000 \
        -c 10 \
        -nt 10 \
        -log F
done

# Lower canine without hominins
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10_no_hominin.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10_no_hominin.tsv" \
    -n "$nc" \
    -p 1000 \
    -s 1000 \
    -c 10 \
    -nt 10 \
    -log F

# Lower canine species means
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10_species_means.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10_species_means.tsv" \
    -n "$nc" \
    -p 1000 \
    -s 1000 \
    -c 10 \
    -nt 10 \
    -log F

# Upper second incisor
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_dec3_10_no_pongo.tsv" \
    -n "$nc" \
    -p 1000 \
    -s 1000 \
    -c 10 \
    -nt 10 \
    -log F
    
# ui2 species means
BURL/build/burl \
    -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo_species_means.csv" \
    -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
    -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_dec3_10_species_means.tsv" \
    -n "$nc" \
    -p 1000 \
    -s 1000 \
    -c 10 \
    -nt 10 \
    -log F

#Gelman-Rubin
for i in 1 2 3 4; do
    BURL/build/burl \
        -i "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo.csv" \
        -it "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt" \
        -o "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/gelmanRubin/out${i}.tsv" \
        -n "$nc" \
        -p 1000 \
        -s 1000 \
        -c 10 \
        -nt 10 \
        -log F
done

# posterior analysis
bash postAnalysis.sh
