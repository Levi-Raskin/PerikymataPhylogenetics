# BURL - Between and within-group Uncertainty in Rates across Lineages.
BURL is a hierarchical Bayesian tool for phylogenetic comparative analysis of multivariate traits. It estimates a unique intraspecific distribution for each taxon, propagating uncertainty in 1) taxon means, 2) taxon variance-covariance matrices, and 3) missing data within a taxon into evolutionary parameter inference. BURL's output can be used in downstream comparative analyses on taxon means, taxon intraspecific covariance structure, or evolutionary rates and covariances, effectively propagating uncertainty in each parameter of interest into a downstream analysis, like phylogenetic ANOVA or modularity analyses.

BURL is described in Raskin et al. (2026) _A hierarchical Bayesian framework accommodates intraspecific and interspecific variation in multivariate traits_. This repo contains all the code needed to reproduce our analyses in that paper, as well as all the command line tools needed to run BURL (our implementation of the hierarchical model described in our paper) on other datasets. Please reach out to [Levi Raskin](mailto:levi_raskin@berkeley.edu) if there are any issues and I would be very happy to help. The data that supports the findings of this study are openly available on [Dryad](https://doi.org/10.5061/dryad.nzs7h4565).
 
To clone the repo and build BURL as a command line tool, run the following bash from GitHub/PerikymataPhylogenetics/:
```{bash}
git clone https://github.com/Levi-Raskin/PerikymataPhylogenetics
cd PerikymataPhylogenetics/scripts/BURL
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

To run an example analysis from command line. Replace filepaths with full-length file paths on your own machine.
```{bash}
./burl \
  -i PerikymataPhylogenetics/data/LCdec3_10.csv \
  -it PerikymataPhylogenetics/data/tree.txt \
  -o PerikymataPhylogenetics/out.tsv \
  -n 1000000 \
  -p 1000 \
  -s 1000 \
  -c 4 \
  -nt 4
```

BURL flags:
* Required
  * -i Input filepath 
  * -it Input tree filepath
  * -o Output filepath
  * -n Chain length; how many cycles do you want to run MCMC for? 
  * -p How often do you want to print useful diagnostics to command line?
  * -s How often do you want to save the current MCMC state as a sample from your posterior distribution
  * -c How many MCMC chains do you want to run? 1 for MCMC, 2+ for Metropolis-coupled MCMC
  * -nt How many threads do you want to use? -nt must be <= -c and at most 1 less than your system's maximum available threads (if you do not know this, BURL will automatically detect maximum threads)
* Optional
  * -h prints help statement and exits program; standalone
  * -log Do you want to log-transform your dataset? This is relatively untested and would make interpreting rates and covariances difficult. This option is included only to facilitate analyses of exceptionally variant datasets, but default is turned off and we recommend it stays off.

Run the following from command line for help:
```{bash}
./burl -h
```

To run BURL on your own data, you need:
* A dataset of multiple continuous traits stored as a .TSV or .CSV file (auto detected). The first row should be a header row with trait names and the first column should be rownames, indicating which taxon a given observation came from. Each row should correspond to a different individual (i.e., row 2 is all the trait values for individual 1). At least one taxon should have multiple individuals sampled.
  * Missing data can be encoded using any non-numeric operator; "NA" is what we used here
  * Missing data is only supported at the moment for taxa where multiple individuals are sampled
  * If your dataset includes only one observation per taxon, BURL is equivalent to standard multivariate Brownian motion with no estimation of intraspecific variation
* A phylogenetic tree in newick form stored as in one line a .txt or .nwk file

This repository has three subfolders:

## data/
* Contains the raw datasets needed to reproduce these analyses, also available on [Dryad](https://doi.org/10.5061/dryad.nzs7h4565).
* pkSpacingCombined.csv is the combined perikymata per millimeter per decile for both the lower canine and upper second incisor
* The other files (e.g., LCdec3_10.csv, UI2dec3_10_no_pongo.csv) are all BURL-friendly formats of the same data that is in pkSpacingCombined.csv
* LCdec3_10_species_means.csv and UI2dec3_10_no_pongo_species_means.csv both encode the empirical species mean perikymata spacing for each taxon in a BURL-friendly format

## figures/
* contains all figure files that appear in the manuscript and supplemental information

## scripts/
* Analysis.R
  * Contains all R-based post processing that's too computationally inexpensive to justify its own command line implementation
  * Effective sample size, KL divergence , symmetrized KL divergence, phylopars, Fréchet variance, modularity analysis, and posterior predictive checks are all here
* Figures.R
  * Contains all the R code we used for figures, except for the figures made in Adobe Illustrator
* fullAnalysis.sh
  * runs all inference and posterior analysis scripts we used to generate our results in our paper
* postAnalysis.sh
  * runs all the post-processing scripts on the posterior distribution
* BURL/
  * C++ source code for BURL
* BURL_coverage/
  * C++ source code to check the coverage of BURL; simulates data under the model and then checks whether the simulated parameter is in the credible interval
* cmd_line_rscripts/
  * R scripts designed to run from the command line to post-process BURL output. See postAnalysis.sh for how to run from command line
  * data_wrangling.R
    * extracts VCV matrices as R matrices from raw BURL posterior output
  * mcmc_convergence.R
    * reads in multiple posterior outputs of the format "filename{1,2,3,4}.tsv" and calculates the rank-normalized Gelman-Rubin statistic to assess MCMC convergence
  * posterior_fits.R
    * Fits an inverse Wishart distribution to the posterior distributions for each VCV matrix using method of moments
  * posterior_pred_draws.R; posterior_pred_draws_ui2.R
    * Draws from the posterior predictive distribution
