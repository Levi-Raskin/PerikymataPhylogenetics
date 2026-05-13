library(coda)
library(data.table)
library(dplyr)
library(MCMCpack)
library(overlapping)
library(parallel)
library(RiemBase)


# functions ---------------------------------------------------------------

calcKLDivergenceInverseWishart <- function(scalePost, dofPost, scalePrior, dofPrior){
  p<- 8
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
calcSymmetrizedKLDivergence <- function(posteriorFit1, posteriorFit2){
  klforward <- calcKLDivergenceInverseWishart(
    scalePost = posteriorFit1$scale,
    dofPost = posteriorFit1$nu,
    scalePrior = posteriorFit2$scale,
    dofPrior = posteriorFit2$nu
  )
  klbackward <- calcKLDivergenceInverseWishart(
    scalePost = posteriorFit2$scale,
    dofPost = posteriorFit2$nu,
    scalePrior = posteriorFit1$scale,
    dofPrior = posteriorFit1$nu
  )
  
  return(
    round(klforward + klbackward, 2)
  )
}



# read data ----------------------------------------------------------------
input <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/"

lc_posterior <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10.tsv"))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

lc_posterior_no_hominin <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_no_hominin.tsv")))
lc_posterior_no_hominin  <- lc_posterior_no_hominin[round(0.1 * nrow(lc_posterior_no_hominin)) : nrow(lc_posterior_no_hominin), ] #apply burnin

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ] #apply burnin


# ESS/GR
lc_gr <- readRDS(paste0(input, "lc/lc_dec3_10_ess_gelman_rubin.RDS"))
ui2_gr <- readRDS(paste0(input, "ui2/ui2_dec3_10_ess_gelman_rubin.RDS"))

#posterior fits
lc_posteriorFits <- readRDS(paste0(input, "lc/lc_dec3_10_posterior_fits.RDS"))
lc_no_hominin_posteriorFits <- readRDS(paste0(input, "lc/lc_dec3_10_no_hominin_posterior_fits.RDS"))
ui2_posteriorFits <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_posterior_fits.RDS"))
ui2_species_means_posteriorFits <- readRDS(paste0(input, "ui2/ui2_dec3_10_posterior_fits_species_means.RDS"))

#posterior predictive
hs_preds <- readRDS(paste0(input, "lc/posteriorPredictive/hsPostPred.rds"))
ne_preds <- readRDS(paste0(input, "lc/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds <- readRDS(paste0(input, "lc/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds <- readRDS(paste0(input, "lc/posteriorPredictive/pantroglodytesPostPred.rds"))
gb_preds <- readRDS(paste0(input, "lc/posteriorPredictive/gorrillaberingeiPostPred.rds"))
gg_preds <- readRDS(paste0(input, "lc/posteriorPredictive/gorillagorillaPostPred.rds"))
pa_preds <- readRDS(paste0(input, "lc/posteriorPredictive/pongoabeliiPostPred.rds"))
ppyg_preds <- readRDS(paste0(input, "lc/posteriorPredictive/pongopygmaeusPostPred.rds"))
hs_preds_ui2 <- readRDS(paste0(input, "ui2/posteriorPredictive/hsPostPred.rds"))
ne_preds_ui2 <- readRDS(paste0(input, "ui2/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds_ui2 <- readRDS(paste0(input, "ui2/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds_ui2 <- readRDS(paste0(input, "ui2/posteriorPredictive/pantroglodytesPostPred.rds"))

#phylopars results
phyloparsRes <-readRDS(paste0(input, "lc/phylopars/lc_dec3_10_phylopars.rds"))

#VCV lists
lc_vcv_list <- readRDS(paste0(input, "lc/lc_dec3_10_vcv_extracted.RDS"))
lc_vcv_list_no_hominins <- readRDS(paste0(input, "lc/lc_dec3_10_no_hominin_vcv_extracted.RDS"))
ui2_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_dec3_10_no_pongo_vcv_extracted.RDS")
lc_vcv_list_species_means <- lc_vcv_list <- readRDS(paste0(input, "lc/lc_dec3_10_vcv_extracted_species_means.RDS"))
ui2_vcv_list_species_means <- lc_vcv_list <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_vcv_extracted_species_means.RDS"))

#Frechet variances
lc_Frechet_Var <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_vcv_frechet_var.RDS")
lc_Frechet_Var_nh <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_vcv_frechet_var_no_hominin.RDS")
ui2_Frechet_Var <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_vcv_frechet_var.RDS")

# ESS/GR ------------------------------------------------------------------
### ess LC
mcmcObj <- mcmc(lc_posterior[,2:ncol(lc_posterior)]) #removes n
ess <- effectiveSize(mcmcObj)
print(ess)
summary(ess)

#GR
summary(unlist(lc_gr))

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


# KL divergence -------------------------------------------------
convertLatexTable <- function(kl, postFits){
  
  order_vec <- c(
    "evolutionary",
    "Homo_sapiens",
    "Neanderthal",
    "Pan_paniscus",
    "Pan_troglodytes",
    "Gorilla_beringei",
    "Gorilla_gorilla",
    "Pongo_abelii",
    "Pongo_pygmaeus"
  )
  
  idx <- match(order_vec, names(kl))
  
  string <- ""
  
  for(i in seq_along(order_vec)){
    
    k <- idx[i]

    kl_val <- if(!is.na(k)) round(kl[k], 2) else "-"
    
    pf_val <- "-"
    if(!is.na(k) && !is.null(postFits[[k]]) && !is.null(postFits[[k]]$nu)){
      pf_val <- round(postFits[[k]]$nu, 2)
    }
    
    string <- paste0(
      string,
      kl_val, " (", pf_val, ")"
    )
    
    if(i != length(order_vec)){
      string <- paste0(string, " & ")
    }
  }
  
  print(order_vec)
  print(string)
}
p <- 8

#### Lower Canine ####

### calc KL divergence
priorDOF<- 10 # numtraits + 2; E[IW] = scale / (dof - p -1 )
priorScale <- matrix(1e-6, p, p)
diag(priorScale) <- 1.0

res <- c()
for(i in 1:length(lc_posteriorFits)){
  res <- c(res, calcKLDivergenceInverseWishart(
    scalePost = lc_posteriorFits[[i]]$scale,
    dofPost = lc_posteriorFits[[i]]$nu,
    scalePrior = priorScale,
    dofPrior = priorDOF
  ))
}
names(res) <- names(lc_posteriorFits)
print(res)
convertLatexTable(res, lc_posteriorFits)

#### Lower Canine no hominins ####
priorDOF<- 10 # numtraits + 2; E[IW] = scale / (dof - p -1 )
priorScale <- matrix(1e-6, p, p)
diag(priorScale) <- 1.0

res <- c()
for(i in 1:length(lc_no_hominin_posteriorFits)){
  res <- c(res, calcKLDivergenceInverseWishart(
    scalePost = lc_no_hominin_posteriorFits[[i]]$scale,
    dofPost = lc_no_hominin_posteriorFits[[i]]$nu,
    scalePrior = priorScale,
    dofPrior = priorDOF
  ))
}
names(res) <- names(lc_no_hominin_posteriorFits)
print(res)
convertLatexTable(res, lc_no_hominin_posteriorFits)

#### Upper second incisor ####
  
### calc KL divergence
priorDOF<- 10 # numtraits + 2; E[IW] = scale / (dof - p -1 )
priorScale <- matrix(1e-6, p, p)
diag(priorScale) <- 1.0

res <- c()
for(i in 1:length(ui2_posteriorFits)){
  res <- c(res, calcKLDivergenceInverseWishart(
    scalePost = ui2_posteriorFits[[i]]$scale,
    dofPost = ui2_posteriorFits[[i]]$nu,
    scalePrior = priorScale,
    dofPrior = priorDOF
  ))
}
names(res) <- names(ui2_posteriorFits)
print(res)
convertLatexTable(res, ui2_posteriorFits)

#### symmetrized KL divergence ####

calcSymmetrizedKLDivergence(lc_posteriorFits$Homo_sapiens, lc_posteriorFits$Neanderthal)
calcSymmetrizedKLDivergence(lc_posteriorFits$Pan_troglodytes, lc_posteriorFits$Pan_paniscus)
calcSymmetrizedKLDivergence(lc_posteriorFits$Gorilla_beringei, lc_posteriorFits$Gorilla_gorilla)
calcSymmetrizedKLDivergence(lc_posteriorFits$Pongo_abelii, lc_posteriorFits$Pongo_pygmaeus)

calcSymmetrizedKLDivergence(ui2_posteriorFits$Homo_sapiens, ui2_posteriorFits$Neanderthal)
calcSymmetrizedKLDivergence(ui2_posteriorFits$Pan_troglodytes, ui2_posteriorFits$Pan_paniscus)

# Rphylopars --------------------------------------------------------------

library(Rphylopars)
dat <- read.csv("Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv")
cn <- colnames(dat)
cn[1] <- "species"
colnames(dat) <- cn

tree <- ape::read.tree(file = "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt")
res <- phylopars(dat, tree)
saveRDS(res, paste0(input, "lc/phylopars/lc_dec3_10_phylopars.rds"))

phyloparsEvoVCV <- phyloparsRes$pars$phylocov
phyloparsIntraVCV <- phyloparsRes$pars$phenocov

#calculate Affine-Invariant Riemannian Metric (AIRM) distnace
library(RiemBase)

nearest_spd <- function(mat, tol = 1e-12) {
  mat <- (mat + t(mat)) / 2
  eig <- eigen(mat, symmetric = TRUE)
  eig$values[eig$values <= 0] <- tol
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

calcAIRMEvoVCV <- function(mat){
  data_list <- riemfactory(list(mat, phyloparsEvoVCV), name = "spd")
  return(rbase.pdist(data_list)[1,2])
}
calcAIRMIntraVCV <- function(mat){
  data_list <- riemfactory(list(mat, phyloparsIntraVCV), name = "spd")
  return(rbase.pdist(data_list)[1,2])
}

for(i in 1:length(lc_vcv_list)){
  mat <- lc_vcv_list[[i]]
  name <- names(lc_vcv_list)[i]
  if(name=="evolutionary"){
    res <- mclapply(
      mat,
      calcAIRMEvoVCV,
      mc.cores = detectCores()-1
    )
  }else{
    res <- mclapply(
      mat,
      calcAIRMIntraVCV,
      mc.cores = detectCores()-1
    )
  }
  saveRDS(res,
          paste0(paste0(input, "lc/phylopars/",
                 name,
                 "_AIRM_distances.rds")))
}

# Exhaustive search through every modularity hypothesis -------------------
lc_evolutionary <- lc_vcv_list$evolutionary
lc_evolutionary_no_hominin <- lc_vcv_list_no_hominins$evolutionary
ui2_evolutionary <- ui2_vcv_list$evolutionary

library(evolqg)
library(parallel)


#### helper functions ####
make_hypo_mat <- function(clustering, n = 8) {
  mat <- matrix(0, n, n)
  for (cluster in clustering) {
    mat[cluster, cluster] <- 1
  }
  mat
}
enumerate_sequential_clusterings <- function(n, K) {
  gaps <- 1:(n - 1)
  cut_combos <- combn(gaps, K - 1, simplify = FALSE)
  lapply(cut_combos, function(cuts) {
    breaks <- c(0, cuts, n)
    lapply(seq_len(K), function(k) (breaks[k] + 1):breaks[k + 1])
  })
}
clustering_label <- function(clustering) {
  paste(sapply(clustering, function(g) paste0(min(g), "-", max(g))), collapse = " | ")
}


#### run analysis ####
n_traits <- 8
all_clusterings <- unlist(
  lapply(2:7, function(K) enumerate_sequential_clusterings(n_traits, K)),
  recursive = FALSE
)

hypo_mats  <- lapply(all_clusterings, make_hypo_mat, n = n_traits)
hypo_labels <- sapply(all_clusterings, clustering_label)
hypo_K      <- sapply(all_clusterings, length)

cat(sprintf("Total hypotheses to test: %d\n", length(hypo_mats)))

calcAVG_ratio <- function(vcvList, permMat) {
  results <- mclapply(vcvList,
                      function(vcv) {
                        cor_mat <- cov2cor(vcv)
                        r <- CalcAVG(permMat, cor_mat)
                        r[1] / r[2]
                      },
                      mc.cores = detectCores() - 1
  )
  results <- unlist(results)
  list(
    mean  = mean(results, trim = 0.005, na.rm = TRUE),
    lower = quantile(results, 0.025,   na.rm = TRUE),
    upper = quantile(results, 0.975,   na.rm = TRUE)
  )
}

run_exhaustive_search <- function(vcvList, label = "dataset") {
  cat(sprintf("\n=== Exhaustive sequential clustering search: %s ===\n", label))
  
  results <- vector("list", length(hypo_mats))
  
  for (i in seq_along(hypo_mats)) {
    r <- calcAVG_ratio(vcvList, hypo_mats[[i]])
    results[[i]] <- data.frame(
      hypothesis = hypo_labels[i],
      K          = hypo_K[i],
      mean_avg   = r$mean,
      lower_95   = r$lower,
      upper_95   = r$upper,
      stringsAsFactors = FALSE
    )
    if (i %% 10 == 0) cat(sprintf("  ... tested %d / %d\n", i, length(hypo_mats)))
  }
  
  df <- do.call(rbind, results)
  df <- df[order(-df$mean_avg), ]
  rownames(df) <- NULL
  df
}

results_hominins    <- run_exhaustive_search(lc_evolutionary,           "with hominins")
results_no_hominins <- run_exhaustive_search(lc_evolutionary_no_hominin, "no hominins")
results_ui2 <- run_exhaustive_search(ui2_evolutionary, "ui2")

print_top <- function(df, n = 10, label = "") {
  cat(sprintf("\nTop %d hypotheses — %s\n", n, label))
  cat(sprintf("%-25s  %4s  %6s  %13s\n", "clustering", "K", "mean", "95% CI"))
  cat(strrep("-", 60), "\n")
  for (i in seq_len(min(n, nrow(df)))) {
    cat(sprintf("%-25s  %4d  %6.3f  (%5.3f, %5.3f)\n",
                df$hypothesis[i], df$K[i], df$mean_avg[i],
                df$lower_95[i], df$upper_95[i]
    ))
  }
}

print_top(results_hominins,    n = 10, label = "with hominins")
print_top(results_no_hominins, n = 10, label = "no hominins")
print_top(results_ui2, n = 10, label = "ui2")

write.csv(results_hominins,    "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/exhaustive_avg_hominins.csv",    row.names = FALSE)
write.csv(results_no_hominins, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/exhaustive_avg_no_hominins.csv", row.names = FALSE)
write.csv(results_ui2, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/exhaustive_avg_ui2.csv", row.names = FALSE)


# comparison against species means inference ----------------------------------------
calcSymmetrizedKLDivergence(lc_posteriorFits$evolutionary, lc_species_means_posteriorFits$evolutionary)
calcSymmetrizedKLDivergence(ui2_species_means_posteriorFits$evolutionary, ui2_posteriorFits$evolutionary)

#### overlap % ####
library(overlapping)

lc_posterior <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10.tsv"))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)):nrow(lc_posterior), ]
lc_posterior_species_means <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10_species_means.tsv"))
lc_posterior_species_means <- lc_posterior_species_means[round(0.1 * nrow(lc_posterior_species_means)):nrow(lc_posterior_species_means), ]

vcv_cols_lc <- grep("^evo_vcv_", colnames(lc_posterior), value = TRUE)
vcv_cols_sm <- grep("^evo_vcv_", colnames(lc_posterior_species_means), value = TRUE)
shared_cols <- intersect(vcv_cols_lc, vcv_cols_sm)

ov_vec <- setNames(numeric(length(shared_cols)), shared_cols)

for (col in shared_cols) {
  ov_result <- overlapping::overlap(
    list(
      lc_posterior[[col]],
      lc_posterior_species_means[[col]]
    )
  )
  ov_vec[col] <- ov_result$OV[[1]]
}

idx <- regmatches(shared_cols, regexpr("\\(\\d+,\\d+\\)", shared_cols))
idx_mat <- do.call(rbind, strsplit(gsub("[()]", "", idx), ","))
rows <- as.integer(idx_mat[, 1])
cols <- as.integer(idx_mat[, 2])
dim_size <- max(rows, cols) + 1  # 0-indexed

ov_matrix <- matrix(NA, nrow = dim_size, ncol = dim_size)
for (k in seq_along(shared_cols)) {
  ov_matrix[rows[k] + 1, cols[k] + 1] <- ov_vec[shared_cols[k]]
}

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ] #apply burnin
ui2_posterior_species_means <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_dec3_10_species_means.tsv"))
ui2_posterior_species_means <- ui2_posterior_species_means[round(0.1 * nrow(ui2_posterior_species_means)):nrow(ui2_posterior_species_means), ]

vcv_cols_ui2 <- grep("^evo_vcv_", colnames(ui2_posterior), value = TRUE)
vcv_cols_sm_ui2 <- grep("^evo_vcv_", colnames(ui2_posterior_species_means), value = TRUE)
shared_cols <- intersect(vcv_cols_lc, vcv_cols_sm)

ov_vec <- setNames(numeric(length(shared_cols)), shared_cols)

for (col in shared_cols) {
  ov_result <- overlapping::overlap(
    list(
      ui2_posterior[[col]],
      ui2_posterior_species_means[[col]]
    )
  )
  ov_vec[col] <- ov_result$OV[[1]]
}

idx <- regmatches(shared_cols, regexpr("\\(\\d+,\\d+\\)", shared_cols))
idx_mat <- do.call(rbind, strsplit(gsub("[()]", "", idx), ","))
rows <- as.integer(idx_mat[, 1])
cols <- as.integer(idx_mat[, 2])
dim_size <- max(rows, cols) + 1  # 0-indexed

ov_matrix_ui2 <- matrix(NA, nrow = dim_size, ncol = dim_size)
for (k in seq_along(shared_cols)) {
  ov_matrix_ui2[rows[k] + 1, cols[k] + 1] <- ov_vec[shared_cols[k]]
}
ov_matrix_ui2

matrix_to_latex <- function(mat) {
  rows <- apply(mat, 1, function(row) {
    paste0("& ", paste(round(row, 2), collapse = " & "), " \\\\")
  })
  cat(paste(rows, collapse = "\n"))
}

matrix_to_latex(100*ov_matrix)
matrix_to_latex(100*ov_matrix_ui2)

#### mean SD diff ####
sd_vec <- setNames(numeric(length(shared_cols)), shared_cols)

for (col in shared_cols) {
  sd_vec[col] <- var(lc_posterior[[col]]) - var(lc_posterior_species_means[[col]])
}

idx <- regmatches(shared_cols, regexpr("\\(\\d+,\\d+\\)", shared_cols))
idx_mat <- do.call(rbind, strsplit(gsub("[()]", "", idx), ","))
rows <- as.integer(idx_mat[, 1])
cols <- as.integer(idx_mat[, 2])
dim_size <- max(rows, cols) + 1  # 0-indexed

sd_matrix <- matrix(NA, nrow = dim_size, ncol = dim_size)
for (k in seq_along(shared_cols)) {
  sd_matrix[rows[k] + 1, cols[k] + 1] <- sd_vec[shared_cols[k]]
}

sd_vec_ui2 <- setNames(numeric(length(shared_cols)), shared_cols)

for (col in shared_cols) {
  sd_vec_ui2[col] <- var(ui2_posterior[[col]]) - var(ui2_posterior_species_means[[col]])
}

idx <- regmatches(shared_cols, regexpr("\\(\\d+,\\d+\\)", shared_cols))
idx_mat <- do.call(rbind, strsplit(gsub("[()]", "", idx), ","))
rows <- as.integer(idx_mat[, 1])
cols <- as.integer(idx_mat[, 2])
dim_size <- max(rows, cols) + 1  # 0-indexed

sd_matrix_ui2 <- matrix(NA, nrow = dim_size, ncol = dim_size)
for (k in seq_along(shared_cols)) {
  sd_matrix_ui2[rows[k] + 1, cols[k] + 1] <- sd_vec_ui2[shared_cols[k]]
}

mean_matrix <- matrix(NA, nrow = dim_size, ncol = dim_size)
for (k in seq_along(shared_cols)) {
  mean_matrix[rows[k] + 1, cols[k] + 1] <- bayestestR::map_estimate(lc_posterior[[shared_cols[k]]])$MAP_Estimate- bayestestR::map_estimate(lc_posterior_species_means[[shared_cols[k]]])$MAP_Estimate
}

mean_matrix_ui2 <- matrix(NA, nrow = dim_size, ncol = dim_size)
for (k in seq_along(shared_cols)) {
  mean_matrix_ui2[rows[k] + 1, cols[k] + 1] <- bayestestR::map_estimate(ui2_posterior[[shared_cols[k]]])$MAP_Estimate - bayestestR::map_estimate(ui2_posterior_species_means[[shared_cols[k]]])$MAP_Estimate
}


matrix_to_latex_mean_sd <- function(mat_mean, mat_sd) {
  rows <- character(nrow(mat_mean))
  for (i in seq_len(nrow(mat_mean))) {
    cells <- paste0(round(mat_mean[i, ], 2), " (", round(mat_sd[i, ], 2), ")")
    rows[i] <- paste0("& ", paste(cells, collapse = " & "), " \\\\")
  }
  cat(paste(rows, collapse = "\n"))
}
matrix_to_latex(mean_matrix)
matrix_to_latex(mean_matrix_ui2)


### Frechet Variance ####
calculateFrechetVariance <- function(vcvList){
  spd_data <- riemfactory(vcvList, name = "spd")
  frechet_mean <- rbase.mean(spd_data)$x
  geodesic_dists <- mclapply(lc_vcv_list$evolutionary, function(mat) {
    pair <- riemfactory(list(mat, frechet_mean), name = "spd")
    rbase.pdist(pair)[1, 2]
  }, mc.cores = detectCores() - 1)
  geodesic_dists <- unlist(geodesic_dists)
  spread <- mean(geodesic_dists^2)
  return(spread)
}

### LC
lc_Frechet_Var <- vector(length = length(lc_vcv_list))
names(lc_Frechet_Var) <- names(lc_vcv_list)
for(i in 1:length(lc_vcv_list)){
  lc_Frechet_Var[i] <- calculateFrechetVariance(lc_vcv_list[[i]])
}
lc_Frechet_Var
saveRDS(lc_Frechet_Var, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_vcv_frechet_var.RDS")

lc_Frechet_Var_nh <- vector(length = length(lc_vcv_list_no_hominins))
names(lc_Frechet_Var_nh) <- names(lc_vcv_list_no_hominins)
for(i in 1:length(lc_vcv_list_no_hominins)){
  lc_Frechet_Var_nh[i] <- calculateFrechetVariance(lc_vcv_list_no_hominins[[i]])
}
lc_Frechet_Var_nh
saveRDS(lc_Frechet_Var_nh, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_vcv_frechet_var_no_hominin.RDS")

### UI2
ui2_Frechet_Var <- vector(length = length(ui2_vcv_list))
names(ui2_Frechet_Var) <- names(ui2_vcv_list)
for(i in 1:length(ui2_vcv_list)){
  ui2_Frechet_Var[i] <- calculateFrechetVariance(ui2_vcv_list[[i]])
}
ui2_Frechet_Var
saveRDS(ui2_Frechet_Var, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_vcv_frechet_var.RDS")

### species means Frechet
print("LC evo")
calculateFrechetVariance(lc_vcv_list$evolutionary)
calculateFrechetVariance(lc_vcv_list_species_means$evolutionary)
print("UI2 evo")
calculateFrechetVariance(ui2_vcv_list$evolutionary)
calculateFrechetVariance(ui2_vcv_list_species_means$evolutionary)


# Posterior mean uncertainty ----------------------------------------------

library(stringr)
matrix_to_latex <- function(mat) {
  rows <- apply(mat, 1, function(row) {
    paste0("& ", paste(round(row, 2), collapse = " & "), " \\\\")
  })
  cat(paste(rows, collapse = "\n"))
}

calculate_posterior_uncertainty <- function(posterior_df) {
  taxon_map <- c(
    "Homo_sapiens"        = "Modern human",
    "Neanderthal"         = "Neandertal",
    "Pan_paniscus"        = "Pan paniscus",
    "Pan_troglodytes"     = "Pan troglodytes",
    "Gorilla_beringei"    = "G. beringei",
    "Gorilla_gorilla"     = "G. gorilla",
    "Pongo_abelii"        = "P. abelii",
    "Pongo_pygmaeus"      = "P. pygmaeus"
  )
  
  present_taxa <- names(taxon_map)[sapply(names(taxon_map), function(taxon) {
    any(str_detect(names(posterior_df), paste0("^", taxon, "_mean_\\d+$")))
  })]
  taxon_map <- taxon_map[present_taxa]
  
  result <- sapply(names(taxon_map), function(taxon) {
    
    taxon_mean_cols <- names(posterior_df)[str_detect(names(posterior_df),
                                                      paste0("^", taxon, "_mean_\\d+$"))]
    decile_idx <- as.integer(str_extract(taxon_mean_cols, "\\d+$"))
    taxon_mean_cols <- taxon_mean_cols[order(decile_idx)]
    
    mean_samples <- as.matrix(posterior_df[, taxon_mean_cols])
    sum(sapply(posterior_df[, taxon_mean_cols], var))
  })
  
  out <- as.data.frame(t(result))
  colnames(out) <- taxon_map[names(taxon_map)]
  out
}

lc_mean_uncertainty <- calculate_posterior_uncertainty(lc_posterior)
print(lc_mean_uncertainty)
matrix_to_latex(lc_mean_uncertainty)

lc_mean_uncertainty_no_hominin <- calculate_posterior_uncertainty(lc_posterior_no_hominin)
print(lc_mean_uncertainty_no_hominin)
matrix_to_latex(lc_mean_uncertainty_no_hominin)

ui2_mean_uncertainty <- calculate_posterior_uncertainty(ui2_posterior)
print(ui2_mean_uncertainty)
matrix_to_latex(ui2_mean_uncertainty)



# relationship of evo VCV to intra VCV ------------------------------------

for(i in 2:length(lc_posteriorFits)){
  print(names(lc_posteriorFits)[i])
  print(calcSymmetrizedKLDivergence(lc_posteriorFits$evolutionary, lc_posteriorFits[[i]]))
  print("======")
}

for(i in 2:length(lc_no_hominin_posteriorFits)){
  print(names(lc_no_hominin_posteriorFits)[i])
  print(calcSymmetrizedKLDivergence(lc_no_hominin_posteriorFits$evolutionary, lc_no_hominin_posteriorFits[[i]]))
  print("======")
}

for(i in 2:length(ui2_posteriorFits)){
  print(names(ui2_posteriorFits)[i])
  print(calcSymmetrizedKLDivergence(ui2_posteriorFits$evolutionary, ui2_posteriorFits[[i]]))
  print("======")
}
