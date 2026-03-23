library(coda)
library(dplyr)
library(MCMCpack)
library(parallel)

lc_posterior <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8.tsv")
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

ui2_posterior <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/ui2_dec3_8.tsv")
ui2_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin


# ESS/GR ------------------------------------------------------------------
### ess LC
mcmcObj <- mcmc(lc_posterior[,2:ncol(lc_posterior)]) #removes n
ess <- effectiveSize(mcmcObj)
print(ess)
summary(ess)

### Gelman Rubin
folder_path <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/GR_test_chains/"  # adjust to your path
file_prefix <- "out"
file_suffix <- ".tsv"
n_chains <- 3
chain_files <- file.path(folder_path, paste0(file_prefix, 0:(n_chains - 1), file_suffix))
chains <- lapply(chain_files, function(file) {
  if (!file.exists(file)) stop(paste("File not found:", file))
  df <- read.table(file, header = TRUE, sep = "\t")
  df <- df[round(0.1*nrow(df)):nrow(df),] #burnin
  mcmc(df)
})
mcmc_list <- mcmc.list(chains)
cat("Gelman-Rubin Diagnostic (R-hat):\n")
gr <- gelman.diag(mcmc_list, autoburnin = TRUE, multivariate = FALSE)
print(gr)
summary(gr$psrf)

### ess UI2
mcmcObj <- mcmc(ui2_posterior[,2:ncol(ui2_posterior)]) #removes n
ess <- effectiveSize(mcmcObj)
print(ess)
summary(ess)

### Gelman Rubin
folder_path <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/GR_test_chains_longer/"  # adjust to your path
file_prefix <- "out"
file_suffix <- ".tsv"
n_chains <- 4
chain_files <- file.path(folder_path, paste0(file_prefix, 0:(n_chains - 1), file_suffix))
chains <- lapply(chain_files, function(file) {
  if (!file.exists(file)) stop(paste("File not found:", file))
  df <- read.table(file, header = TRUE, sep = "\t")
  df <- df[round(0.1*nrow(df)):nrow(df),] #burnin
  mcmc(df)
})
mcmc_list <- mcmc.list(chains)
cat("Gelman-Rubin Diagnostic (R-hat):\n")
gr <- gelman.diag(mcmc_list, autoburnin = TRUE, multivariate = FALSE)
print(gr)
summary(gr$psrf)

# KL divergence -------------------------------------------------

#data wrangling
# extract_vcv <- function(row, prefix, p = 8) {
#   mat <- matrix(NA, p, p)
#   for (i in 0:(p-1)) {
#     for (j in 0:(p-1)) {
#       col <- paste0(prefix, "vcv_.", i, ".", j, ".")
#       if (col %in% names(row)) {
#         mat[i+1, j+1] <- row[[col]]
#       }
#     }
#   }
#   mat
# }
# 
# species_list <- c(
#   "evolutionary",
#   "Pongo_abelii",
#   "Pongo_pygmaeus",
#   "Pan_troglodytes",
#   "Pan_paniscus",
#   "Gorilla_beringei",
#   "Gorilla_gorilla",
#   "Homo_sapiens",
#   "Neanderthal"
# )
# 
# prefix_map <- c(
#   evolutionary    = "evo_",
#   Pongo_abelii    = "Pongo_abelii_",
#   Pongo_pygmaeus  = "Pongo_pygmaeus_",
#   Pan_troglodytes = "Pan_troglodytes_",
#   Pan_paniscus    = "Pan_paniscus_",
#   Gorilla_beringei = "Gorilla_beringei_",
#   Gorilla_gorilla  = "Gorilla_gorilla_",
#   Homo_sapiens     = "Homo_sapiens_",
#   Neanderthal      = "Neanderthal_"
# )
# 
# n_samples <- nrow(lc_posterior)
# 
# vcv_list <- mclapply(species_list, function(sp) {
#   prefix <- prefix_map[sp]
#   lapply(seq_len(n_samples), function(i) {
#     extract_vcv(lc_posterior[i, ], prefix, p = 8)
#   })
# },
# mc.cores = detectCores() -1)
# 
# names(vcv_list) <- species_list
# 
# saveRDS(vcv_list, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_VCVs_extracted.rds")
vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_VCVs_extracted.rds")

p <- 8
dof <- 10 # numtraits + 2; E[IW] = scale / (dof - p -1 )
scale <- matrix(1e-6, p, p)
diag(scale) <- 1.0

lnl <- function(dat, S, nu){
  res <- mclapply(
    dat,
    FUN = function(x) {
      LaplacesDemon::dinvwishart(x, nu, S, log = TRUE)
    },
    mc.cores = detectCores() - 1
  )
  vals <- unlist(res)
  if (any(!is.finite(vals))) return(-Inf)
  sum(vals)
}

posteriorFits <- list()
for(i in 1:length(vcv_list)){
  post <- vcv_list[[i]]
  
  nu_0 <- dof
  S_0  <- scale
  L_0  <- t(chol(S_0 + diag(1e-6, p)))
  
  init <- c(log(nu_0 - p - 1), L_0[lower.tri(L_0, diag = TRUE)])
  
  neg_lnl <- function(params) {
    nu <- exp(params[1]) + p + 1
    
    L <- matrix(0, p, p)
    L[lower.tri(L, diag = TRUE)] <- params[2:(p*(p+1)/2 + 1)]
    S <- L %*% t(L)
    
    eig <- eigen(S, only.values = TRUE)$values
    if (any(eig <= 1e-10)) return(.Machine$double.xmax)
    
    ll <- lnl(post, S, nu)
    if (!is.finite(ll)) return(.Machine$double.xmax)
    -ll
  }
  
  cat("Initial neg log-lik:", neg_lnl(init), "\n")
  
  tictoc::tic("Optimization")
  fit <- optim(
    par     = init,
    fn      = neg_lnl,
    method  = "BFGS",
    control = list(maxit = 1000, reltol = 1e-8, trace = 1)
  )
  tictoc::toc()
  
  nu_hat <- exp(fit$par[1]) + p + 1
  L_hat  <- matrix(0, p, p)
  L_hat[lower.tri(L_hat, diag = TRUE)] <- fit$par[2:(p*(p+1)/2 + 1)]
  S_hat  <- L_hat %*% t(L_hat)
  
  posteriorFits[[i]] <- list(
    nu = nu_hat,
    scale = S_hat
  )
}
names(posteriorFits) <- names(vcv_list)
saveRDS(posteriorFits, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_posterior_fits.rds")

calcKLDivergence <- function(scalePost, dofPost, scalePrior, dofPrior){
  
}


# Analyses on the posterior predictive distributions ----------------------------------------------

hs_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/hsPostPred.rds")
ne_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/neanderthalPostPred.rds")
pp_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/panpaniscusPostPred.rds")
pt_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pantroglodytesPostPred.rds")
gb_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/gorrillaberingeiPostPred.rds")
gg_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/gorillagorillaPostPred.rds")
pa_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pongoabeliiPostPred.rds")
ppyg_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pongopygmaeusPostPred.rds")

### variance in the posterior predictive
for(i in 1:8){
  print(colnames(hs_preds)[i])
  print("Modern human var: ")
  print(var(hs_preds[,i]))
  print("Neanderthal var: ")
  print(var(ne_preds[,i]))
  print("=======")
}
