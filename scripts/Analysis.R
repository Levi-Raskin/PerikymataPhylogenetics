library(coda)
library(dplyr)
library(MCMCpack)
library(overlapping)
library(parallel)

lc_posterior <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8.tsv")
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

lc_posterior_no_hominin <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8_no_hominin.tsv")
lc_posterior_no_hominin  <- lc_posterior_no_hominin[round(0.1 * nrow(lc_posterior_no_hominin)) : nrow(lc_posterior_no_hominin), ] #apply burnin

ui2_posterior <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/ui2_dec3_8.tsv")
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ] #apply burnin



# Functions ---------------------------------------------------------------
convertLatexTable <- function(vec){
  string <- ""
  for(i in vec){
    string <- paste0(string, round(i, 2), " & ")
  }
  print(string)
}
extract_vcv <- function(row, prefix, p = 8) {
  mat <- matrix(NA, p, p)
  for (i in 0:(p-1)) {
    for (j in 0:(p-1)) {
      col <- paste0(prefix, "vcv_.", i, ".", j, ".")
      if (col %in% names(row)) {
        mat[i+1, j+1] <- row[[col]]
      }
    }
  }
  mat
}
calcKLDivergenceInverseWishart <- function(scalePost, dofPost, scalePrior, dofPrior){
  V1 <- solve(scalePost)
  V2 <- solve(scalePrior)
  n1 <- dofPost
  n2 <- dofPrior
  term1 <- n2 * as.numeric(
    determinant(V2, logarithm = TRUE)$modulus - 
      determinant(V1, logarithm = TRUE)$modulus
  )
  term2 <- n1 * sum(diag(solve(V2) %*% V1))
  term3 <- 2 * (CholWishart::lmvgamma(n2/2, p) - CholWishart::lmvgamma(n1/2, p))
  term4 <- (n1 - n2) * CholWishart::mvdigamma(n1/2, p)
  term5 <- -n1 * p
  
  result <- 0.5 * (term1 + term2 + term3 + term4 + term5)
  return(as.numeric(result))
}

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

### ess LC w/o hominin
mcmcObj <- mcmc(lc_posterior_no_hominin[,2:ncol(lc_posterior_no_hominin)]) #removes n
ess <- effectiveSize(mcmcObj)
print(ess)
summary(ess)


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

#### data wrangling ####
# # LC
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
# 


# #LC no hominins
# species_list <- c(
#   "evolutionary",
#   "Pongo_abelii",
#   "Pongo_pygmaeus",
#   "Pan_troglodytes",
#   "Pan_paniscus",
#   "Gorilla_beringei",
#   "Gorilla_gorilla"
# )
# 
# prefix_map <- c(
#   evolutionary    = "evo_",
#   Pongo_abelii    = "Pongo_abelii_",
#   Pongo_pygmaeus  = "Pongo_pygmaeus_",
#   Pan_troglodytes = "Pan_troglodytes_",
#   Pan_paniscus    = "Pan_paniscus_",
#   Gorilla_beringei = "Gorilla_beringei_",
#   Gorilla_gorilla  = "Gorilla_gorilla_"
# )
# 
# n_samples <- nrow(lc_posterior_no_hominin)
# 
# vcv_list <- mclapply(species_list, function(sp) {
#   prefix <- prefix_map[sp]
#   lapply(seq_len(n_samples), function(i) {
#     extract_vcv(lc_posterior_no_hominin[i, ], prefix, p = 8)
#   })
# },
# mc.cores = detectCores() -1)
# 
# names(vcv_list) <- species_list
# saveRDS(vcv_list, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_no_hominin_VCVs_extracted.rds")

# UI2
# species_list <- c(
#   "evolutionary",
#   "Pan_troglodytes",
#   "Pan_paniscus",
#   "Homo_sapiens",
#   "Neanderthal"
# )
# 
# prefix_map <- c(
#   evolutionary    = "evo_",
#   Pan_troglodytes = "Pan_troglodytes_",
#   Pan_paniscus    = "Pan_paniscus_",
#   Homo_sapiens     = "Homo_sapiens_",
#   Neanderthal      = "Neanderthal_"
# )
# 
# n_samples <- nrow(ui2_posterior)
# 
# vcv_list <- mclapply(species_list, function(sp) {
#   prefix <- prefix_map[sp]
#   lapply(seq_len(n_samples), function(i) {
#     extract_vcv(ui2_posterior[i, ], prefix, p = 8)
#   })
# },
# mc.cores = detectCores() -1)
# 
# names(vcv_list) <- species_list
# 
# saveRDS(vcv_list, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/ui2_VCVs_extracted.rds")

#### Lower Canine ####
lc_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_VCVs_extracted.rds")

p <- 8

### MLE fitting to posterior VCV distributions
posteriorFits <- list()
#moment matching
for (i in 1:length(lc_vcv_list)) {
  post <- lc_vcv_list[[i]]
  N <- length(post)
  
  M <- Reduce("+", lapply(post, solve)) / N

  f <- function(nu) {
      determinant(M, logarithm = T)$modulus - p * log(nu / 2) +
      sum(digamma((nu - 0:(p - 1)) / 2)) -
      mean(sapply(post, function(x) determinant(solve(x), logarithm = T)$modulus))
  }
  
  f_prime <- function(nu) {
    -p / nu + 0.5 * sum(trigamma((nu - 0:(p - 1)) / 2))
  }
  
  nu <- p + 2
  while (f(nu) >= 0) nu <- p + 1 + (nu - p - 1) / 2
  
  for (iter in 1:100) {
    step <- f(nu) / f_prime(nu)
    nu <- nu - step
    if (abs(step) < 1e-12) break
  }
  
  S <- nu * solve(M)
  
  posteriorFits[[i]] <- list(
    nu = as.numeric(nu),
    scale = S
  )
}
names(posteriorFits) <- names(lc_vcv_list)
saveRDS(posteriorFits, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_posterior_fits.rds")
posteriorFits <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_posterior_fits.rds")

### calc KL divergence
priorDOF<- 10 # numtraits + 2; E[IW] = scale / (dof - p -1 )
priorScale <- matrix(1e-6, p, p)
diag(priorScale) <- 1.0

res <- c()
for(i in 1:length(posteriorFits)){
  res <- c(res, calcKLDivergenceInverseWishart(
    scalePost = posteriorFits[[i]]$scale,
    dofPost = posteriorFits[[i]]$nu,
    scalePrior = priorScale,
    dofPrior = priorDOF
  ))
}
names(res) <- names(posteriorFits)
print(res)
convertLatexTable(res)

for(i in 1:length(posteriorFits)){
  print(paste(names(posteriorFits)[i],round( posteriorFits[[i]]$nu , 2)))
}

#### Lower Canine no hominins ####
lc_no_hominin_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_no_hominin_VCVs_extracted.rds")

p <- 8

### MLE fitting to posterior VCV distributions
posteriorFits <- list()
#moment matching
for (i in 1:length(lc_no_hominin_vcv_list)) {
  post <- lc_no_hominin_vcv_list[[i]]
  N <- length(post)
  
  M <- Reduce("+", lapply(post, solve)) / N
  
  f <- function(nu) {
    determinant(M, logarithm = T)$modulus - p * log(nu / 2) +
      sum(digamma((nu - 0:(p - 1)) / 2)) -
      mean(sapply(post, function(x) determinant(solve(x), logarithm = T)$modulus))
  }
  
  f_prime <- function(nu) {
    -p / nu + 0.5 * sum(trigamma((nu - 0:(p - 1)) / 2))
  }
  
  nu <- p + 2
  while (f(nu) >= 0) nu <- p + 1 + (nu - p - 1) / 2
  
  for (iter in 1:100) {
    step <- f(nu) / f_prime(nu)
    nu <- nu - step
    if (abs(step) < 1e-12) break
  }
  
  S <- nu * solve(M)
  
  posteriorFits[[i]] <- list(
    nu = as.numeric(nu),
    scale = S
  )
}
names(posteriorFits) <- names(lc_no_hominin_vcv_list)
saveRDS(posteriorFits, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_no_hominin_posterior_fits.rds")
posteriorFits <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_no_hominin_posterior_fits.rds")

### calc KL divergence
priorDOF<- 10 # numtraits + 2; E[IW] = scale / (dof - p -1 )
priorScale <- matrix(1e-6, p, p)
diag(priorScale) <- 1.0

res <- c()
for(i in 1:length(posteriorFits)){
  res <- c(res, calcKLDivergenceInverseWishart(
    scalePost = posteriorFits[[i]]$scale,
    dofPost = posteriorFits[[i]]$nu,
    scalePrior = priorScale,
    dofPrior = priorDOF
  ))
}
names(res) <- names(posteriorFits)
print(res)
convertLatexTable(res)

for(i in 1:length(posteriorFits)){
  print(paste(names(posteriorFits)[i],round( posteriorFits[[i]]$nu , 2)))
}

#### Upper second incisor ####
ui2_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/ui2_VCVs_extracted.rds")

### MLE fitting to posterior VCV distributions
posteriorFits <- list()

for (i in 1:length(ui2_vcv_list)) {
  post <- ui2_vcv_list[[i]]
  N <- length(post)
  
  P <- lapply(post, solve)
  m_P <- Reduce("+", P) / N
  m_log_det_P <- mean(sapply(P, function(x) as.numeric(determinant(x, logarithm = T)$modulus)))
  log_det_m_P <- as.numeric(determinant(m_P, logarithm = T)$modulus)
  
  f <- function(nu) {
    p * log(2) + log_det_m_P - p * log(nu) +
      sum(digamma((nu + 1 - 1:p) / 2)) - m_log_det_P
  }
  
  f_prime <- function(nu) {
    -p / nu + 0.5 * sum(trigamma((nu + 1 - 1:p) / 2))
  }
  
  nu <- p + 2
  epsilon <- 1e-8
  
  lb <- p - 1 + epsilon
  
  for (iter in 1:100) {
    step <- f(nu) / f_prime(nu)
    nu_s <- nu - step
    if (nu_s <= lb) nu_s <- lb + (nu - lb) / 2
    if (abs(nu_s - nu) < 1e-12) { nu <- nu_s; break }
    nu <- nu_s
  }
  
  V <- nu * solve(m_P)
  
  posteriorFits[[i]] <- list(
    nu = nu, 
    scale = V
  )
}

names(posteriorFits) <- names(ui2_vcv_list)
saveRDS(posteriorFits, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/ui2_posterior_fits.rds")
posteriorFits <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/ui2_posterior_fits.rds")
  
### calc KL divergence
priorDOF<- 10 # numtraits + 2; E[IW] = scale / (dof - p -1 )
priorScale <- matrix(1e-6, p, p)
diag(priorScale) <- 1.0

res <- c()
for(i in 1:length(posteriorFits)){
  res <- c(res, calcKLDivergenceInverseWishart(
    scalePost = posteriorFits[[i]]$scale,
    dofPost = posteriorFits[[i]]$nu,
    scalePrior = priorScale,
    dofPrior = priorDOF
  ))
}
names(res) <- names(posteriorFits)
print(res)
convertLatexTable(res)

for(i in 1:length(posteriorFits)){
  print(paste(names(posteriorFits)[i],round( posteriorFits[[i]]$nu , 2)))
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