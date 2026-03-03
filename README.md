BURL - Between and within-group Uncertainty in Rates across Lineages

Code and repository for Raskin et al. (202X).

To build BURL as a commmand line tool run the following bash from GitHub/PerikymataPhylogenetics/:
```{bash}
cd scripts/BURL
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

And run from command line:
```{bash}
./burl -i /Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv -it /Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt -o /Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/out.tsv -n 100000000 -p 1000 -s 1000 -c 10 -nt 10
```

Run from command line for help:
```{bash}
./burl -h
```
