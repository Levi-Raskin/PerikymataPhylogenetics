library(coda)
library(dplyr)
library(entropy)
library(LaplacesDemon)
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
folder_path <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/GR_test_chains"  # adjust to your path
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

### ess UI2
mcmcObj <- mcmc(ui2_posterior[,2:ncol(ui2_posterior)]) #removes n
ess <- effectiveSize(mcmcObj)
print(ess)
summary(ess)

### Gelman Rubin
folder_path <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/ui2/GR_test_chains"  # adjust to your path
file_prefix <- "chain"
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

# Empirical KL divergence -------------------------------------------------

### tip VCVs
n_traits <- 8
species_list <- c("Pongo_abelii", "Pongo_pygmaeus", "Pan_troglodytes", 
                  "Pan_paniscus", "Gorilla_beringei", "Gorilla_gorilla",
                  "Homo_sapiens", "Neanderthal")

# Draw from NIW prior (tips & evo vcv share same prior)
draw_iw_prior <- function(n_draws, n_traits = 8) {
  nu     <- n_traits + 2
  lambda <- nu - n_traits - 1
  mu0    <- rep(0, n_traits)
  Psi    <- matrix(1e-6, n_traits, n_traits)
  diag(Psi) <- 1.0
  
  mclapply(1:n_draws, function(s) {
    Sigma <- LaplacesDemon::rinvwishart(nu, Psi)
    list(Sigma = Sigma)
  }, mc.cores = detectCores()-1)
}

prior_draws <- draw_iw_prior(n_draws = nrow(lc_posterior))

prior_vcv_mat <- do.call(rbind, lapply(prior_draws, function(d) as.vector(d$Sigma)))

kl_vcv_raw <- mclapply(species_list, function(sp) {
  vcv_cols <- as.vector(outer(0:(n_traits-1), 0:(n_traits-1),
                              function(i, j) paste0(sp, "_vcv_.", i, ".", j, ".")))
  post_vcv <- as.matrix(dplyr::select(lc_posterior, all_of(vcv_cols)))
  
  kl_per_dim <- sapply(1:ncol(post_vcv), function(d) {
    breaks <- seq(
      min(c(post_vcv[, d], prior_vcv_mat[, d])),
      max(c(post_vcv[, d], prior_vcv_mat[, d])),
      length.out = 100
    )
    p <- hist(post_vcv[, d],      breaks = breaks, plot = FALSE)$counts + 1
    q <- hist(prior_vcv_mat[, d], breaks = breaks, plot = FALSE)$counts + 1
    KL.empirical(p, q)
  })
  
  per_dim_df <- data.frame(
    species = sp,
    row     = rep(1:n_traits, each  = n_traits),
    col     = rep(1:n_traits, times = n_traits),
    element = vcv_cols,
    kl      = kl_per_dim
  )
  
  summary_df <- data.frame(
    species     = sp,
    kl_vcv      = sum(kl_per_dim),
    kl_vcv_mean = mean(kl_per_dim)
  )
  
  list(summary = summary_df, per_dim = per_dim_df)
}, mc.cores = detectCores()-1)

kl_vcv_summary <- lapply(kl_vcv_raw, `[[`, "summary") |> bind_rows()
kl_vcv_per_dim <- lapply(kl_vcv_raw, `[[`, "per_dim") |> bind_rows()

kl_vcv_summary
kl_vcv_per_dim

saveRDS(kl_vcv_raw, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_tipVCV_KL_div.rds")

### evolutionary VCV
evo_vcv_cols <- as.vector(outer(0:(n_traits-1), 0:(n_traits-1),
                                function(i, j) paste0("evo_vcv_.", i, ".", j, ".")))

post_evo_vcv <- as.matrix(dplyr::select(lc_posterior, all_of(evo_vcv_cols)))

kl_per_dim <- sapply(1:ncol(post_evo_vcv), function(d) {
  breaks <- seq(
    min(c(post_evo_vcv[, d], prior_vcv_mat[, d])),
    max(c(post_evo_vcv[, d], prior_vcv_mat[, d])),
    length.out = 100
  )
  
  p <- hist(post_evo_vcv[, d],  breaks = breaks, plot = FALSE)$counts + 1  # +1 Laplace smoothing
  q <- hist(prior_vcv_mat[, d], breaks = breaks, plot = FALSE)$counts + 1
  
  KL.empirical(p, q)
})

kl_evo_vcv <- data.frame(
  element     = evo_vcv_cols,
  row         = rep(1:n_traits, each = n_traits),
  col         = rep(1:n_traits, times = n_traits),
  kl          = kl_per_dim
)

kl_evo_vcv
saveRDS(kl_evo_vcv, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_evoVCV_KL_div.rds")

