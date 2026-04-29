library(data.table)
library(dplyr)
library(MASS)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

get_posterior_predictive <- function(posterior, species, n_traits, n_samples) {
  
  mean_cols <- paste0(species, "_mean_", 0:(n_traits - 1))
  mu_samples <- as.matrix(posterior[, mean_cols])
  
  vcv_cols <- outer(0:(n_traits - 1), 0:(n_traits - 1),
                    FUN = function(i, j) paste0(species, "_vcv_(", i, ",", j, ")"))
  vcv_mat <- as.matrix(posterior[, as.vector(vcv_cols)])
  
  draw_one <- function(s) {
    Sigma <- matrix(vcv_mat[s, ], nrow = n_traits, ncol = n_traits)
    Sigma <- (Sigma + t(Sigma)) / 2
    Sigma <- Sigma + diag(1e-6, n_traits)
    mvrnorm(n = 1, mu = mu_samples[s, ], Sigma = Sigma)
  }
  
  preds <- do.call(rbind, mclapply(1:n_samples, draw_one, mc.cores = detectCores() - 1))
  
  colnames(preds) <- trait_labels
  as.data.frame(preds) |>
    mutate(species = species)
}

posterior <- as.data.frame(fread(input))
posterior <- posterior[round(0.1 * nrow(posterior)) : nrow(posterior), ] #apply burnin

n_samples <- nrow(posterior)
n_traits <- 8
trait_labels <- paste0("Decile ", 3:10)

# Generate posterior predictives
hs_preds <- get_posterior_predictive(posterior, "Homo_sapiens", n_traits, n_samples)
ne_preds <- get_posterior_predictive(posterior, "Neanderthal", n_traits, n_samples)
pp_preds <- get_posterior_predictive(posterior, "Pan_paniscus", n_traits, n_samples)
pt_preds <- get_posterior_predictive(posterior, "Pan_troglodytes", n_traits, n_samples)
gb_preds <- get_posterior_predictive(posterior, "Gorilla_beringei", n_traits, n_samples)
gg_preds <- get_posterior_predictive(posterior, "Gorilla_gorilla", n_traits, n_samples)
pa_preds <- get_posterior_predictive(posterior, "Pongo_abelii", n_traits, n_samples)
ppyg_preds <- get_posterior_predictive(posterior, "Pongo_pygmaeus", n_traits, n_samples)
saveRDS(hs_preds, paste0(output, "hsPostPred.rds"))
saveRDS(ne_preds, paste0(output, "neanderthalPostPred.rds"))
saveRDS(pp_preds, paste0(output, "panpaniscusPostPred.rds"))
saveRDS(pt_preds, paste0(output, "pantroglodytesPostPred.rds"))
saveRDS(gb_preds, paste0(output, "gorrillaberingeiPostPred.rds"))
saveRDS(gg_preds, paste0(output, "gorillagorillaPostPred.rds"))
saveRDS(pa_preds, paste0(output, "pongoabeliiPostPred.rds"))
saveRDS(ppyg_preds, paste0(output, "pongopygmaeusPostPred.rds"))

message("----------------------------------------")
message(paste0("Posterior pred draws finished for: ", input))
message("----------------------------------------")