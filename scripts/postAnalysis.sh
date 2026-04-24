#!/bin/bash

RESULTS=/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs

# calculate Rank-Normalized R hat
Rscript cmd_line_rscripts/mcmc_convergence.R \
    $RESULTS/lc/gelmanRubin/out \
    $RESULTS/lc/lc_dec3_10_ess_gelman_rubin.RDS

Rscript cmd_line_rscripts/mcmc_convergence.R \
    $RESULTS/ui2/gelmanRubin/out \
    $RESULTS/ui2/ui2_dec3_10_ess_gelman_rubin.RDS

# wrangling data
Rscript cmd_line_rscripts/data_wrangling.R \
    $RESULTS/lc/lc_dec3_10.tsv \
    $RESULTS/lc/lc_dec3_10_vcv_extracted.RDS

Rscript cmd_line_rscripts/data_wrangling.R \
    $RESULTS/lc/lc_dec3_10_no_hominin.tsv \
    $RESULTS/lc/lc_dec3_10_no_hominin_vcv_extracted.RDS

Rscript cmd_line_rscripts/data_wrangling.R \
    $RESULTS/ui2/ui2_dec3_10_no_pongo.tsv \
    $RESULTS/ui2/ui2_dec3_10_no_pongo_vcv_extracted.RDS

# fit IW to VCV posterior distributions
Rscript cmd_line_rscripts/posterior_fits.R \
    $RESULTS/lc/lc_dec3_10_vcv_extracted.RDS \
    $RESULTS/lc/lc_dec3_10_posterior_fits.RDS

Rscript cmd_line_rscripts/posterior_fits.R \
    $RESULTS/lc/lc_dec3_10_no_hominin_vcv_extracted.RDS \
    $RESULTS/lc/lc_dec3_10_no_hominin_posterior_fits.RDS

Rscript cmd_line_rscripts/posterior_fits.R \
    $RESULTS/ui2/ui2_dec3_10_no_pongo_vcv_extracted.RDS \
    $RESULTS/ui2/ui2_dec3_10_no_pongo_posterior_fits.RDS

# posterior predictive draws (just done on lower canine full dataset)
Rscript cmd_line_rscripts/posterior_pred_draws.R \
    $RESULTS/lc/lc_dec3_10.tsv \
    $RESULTS/lc/posteriorPredictive/
