library(coda)
library(data.table)
library(dplyr)
library(MCMCpack)
library(overlapping)
library(parallel)

input <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/"

lc_posterior <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10.tsv"))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

lc_posterior_no_hominin <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_no_hominin.tsv")))
lc_posterior_no_hominin  <- lc_posterior_no_hominin[round(0.1 * nrow(lc_posterior_no_hominin)) : nrow(lc_posterior_no_hominin), ] #apply burnin

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ] #apply burnin

# Functions ---------------------------------------------------------------
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

#GR
lc_gr <- readRDS(paste0(input, "lc/lc_dec3_10_ess_gelman_rubin.RDS"))

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

ui2_gr <- readRDS(paste0(input, "ui2/ui2_dec3_10_ess_gelman_rubin.RDS"))

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
lc_posteriorFits <- readRDS(paste0(input, "lc/lc_dec3_10_posterior_fits.RDS"))

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
lc_no_hominin_posteriorFits <- readRDS(paste0(input, "lc/lc_dec3_10_no_hominin_posterior_fits.RDS"))
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
ui2_posteriorFits <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_posterior_fits.RDS"))
  
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

#### Lower canine symmetrized KL divergence ####
lc_posteriorFits <- readRDS(paste0(input, "lc/lc_dec3_10_posterior_fits.RDS"))

calcSymmetrizedKLDivergence <- function(posteriorFit1, posteriorFit2){
  p <- 8
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

calcSymmetrizedKLDivergence(lc_posteriorFits$Homo_sapiens, lc_posteriorFits$Neanderthal)
calcSymmetrizedKLDivergence(lc_posteriorFits$Pan_troglodytes, lc_posteriorFits$Pan_paniscus)
calcSymmetrizedKLDivergence(lc_posteriorFits$Gorilla_beringei, lc_posteriorFits$Gorilla_gorilla)
calcSymmetrizedKLDivergence(lc_posteriorFits$Pongo_abelii, lc_posteriorFits$Pongo_pygmaeus)

# Analyses on the posterior predictive distributions ----------------------------------------------

hs_preds <- readRDS(paste0(input, "lc/posteriorPredictive/hsPostPred.rds"))
ne_preds <- readRDS(paste0(input, "lc/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds <- readRDS(paste0(input, "lc/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds <- readRDS(paste0(input, "lc/posteriorPredictive/pantroglodytesPostPred.rds"))
gb_preds <- readRDS(paste0(input, "lc/posteriorPredictive/gorrillaberingeiPostPred.rds"))
gg_preds <- readRDS(paste0(input, "lc/posteriorPredictive/gorillagorillaPostPred.rds"))
pa_preds <- readRDS(paste0(input, "lc/posteriorPredictive/pongoabeliiPostPred.rds"))
ppyg_preds <- readRDS(paste0(input, "lc/posteriorPredictive/pongopygmaeusPostPred.rds"))

### variance in the posterior predictive
for(i in 1:8){
  print(colnames(hs_preds)[i])
  print("Modern human var: ")
  print(var(hs_preds[,i]))
  print("Neanderthal var: ")
  print(var(ne_preds[,i]))
  print("=======")
}

# Modularity test ---------------------------------------------------------
lc_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10_vcv_extracted.RDS")
lc_evolutionary <- lc_vcv_list$evolutionary
lc_vcv_list_no_hominins <- readRDS(paste0(input, "lc/lc_dec3_10_no_hominin_vcv_extracted.RDS"))
lc_evolutionary_no_hominin <- lc_vcv_list_no_hominins$evolutionary

#### AVG Ratio #### 
library(evolqg)
testAVG <- function(vcvList, permMat){
  avg_rat_results <- mclapply(vcvList, 
                              function(vcv) {
                                cor_mat <- cov2cor(vcv)
                                result <- CalcAVG(permMat, cor_mat)
                                return(result[1] /result[2])
                              },
                              mc.cores = detectCores()-1)
  
  avg_rat_results <- unlist(avg_rat_results)
  
  meanAvg <-   mean(avg_rat_results, trim = 0.005, na.rm = TRUE)
  quant <-quantile(avg_rat_results, c(0.025, 0.975), na.rm = TRUE)
  
  print(paste0(round(meanAvg, 2), 
              " (", 
               round(quant[1],2), 
               ", ", 
               round(quant[2],2), ")" ))
}

hypoMat1 <- matrix(0, 8, 8)
hypoMat1[1:4, 1:4] <- 1
hypoMat1[5:8, 5:8] <- 1

hypoMat2 <- matrix(0, 8, 8)
hypoMat2[1:3, 1:3] <- 1
hypoMat2[4:8, 4:8] <- 1

hypoMat3 <- matrix(0, 8, 8)
hypoMat3[1:3, 1:3] <- 1
hypoMat3[4,4] <- 1
hypoMat3[5:8, 5:8] <- 1

#### C1 analysis with hominins ####

#### Hypothesis 1: deciles 3-6; deciles 7-10 ####
testAVG(lc_evolutionary, hypoMat1)

#### Hypothesis 2: deciles 3-5; deciles 6-10 ####
testAVG(lc_evolutionary, hypoMat2)

#### Hypothesis 3: deciles 3-5; decile 6; deciles 7-10 ####
testAVG(lc_evolutionary, hypoMat3)

#### C1 analysis without hominins ####
#### Hypothesis 1: deciles 3-6; deciles 7-10 ####
testAVG(lc_evolutionary_no_hominin, hypoMat1)

#### Hypothesis 2: deciles 3-5; deciles 6-10 ####
testAVG(lc_evolutionary_no_hominin, hypoMat2)

#### Hypothesis 3: deciles 3-5; decile 6; deciles 7-10 ####
testAVG(lc_evolutionary_no_hominin, hypoMat3)

#### CR ####
CalcCR <- function(pv, vcv_matrix) {
  ### modified from Adams (2016)
  gps<-factor(pv)
  S11<-S11.0<-vcv_matrix[which(gps==levels(gps)[1]),which(gps==levels(gps)[1])]
  S22<-S22.0<-vcv_matrix[which(gps==levels(gps)[2]),which(gps==levels(gps)[2])]
  diag(S11.0)<-0
  diag(S22.0)<-0
  S12<-vcv_matrix[which(gps==levels(gps)[1]),which(gps==levels(gps)[2])]
  S21<-t(S12)
  return(sqrt( sum(diag(S12%*%S21)) / sqrt(sum(diag(S11.0%*%S11.0))*sum(diag(S22.0%*%S22.0)))))
}
testCR <- function(vcvList, pm){
  cr_results <- mclapply(vcvList, 
                         function(vcv) CalcCR(pm, vcv),
                         mc.cores = detectCores()-1)
  cr_results <- unlist(cr_results)
  cat("Posterior mean CR:", mean(cr_results, na.rm = TRUE), "\n")
  cat("95% credible interval:", quantile(cr_results, c(0.025, 0.975), na.rm = TRUE), "\n")
}

hypo1Vec <- c(rep(1, 4), rep (2, 4))
hypo2Vec <- c(rep(1, 3), rep (2, 5))

#### C1 analysis with hominins ####
#### Hypothesis 1: deciles 3-6; deciles 7-10 ####
testCR(lc_evolutionary, hypo1Vec)

#### Hypothesis 2: deciles 3-5; deciles 6-10 ####
testCR(lc_evolutionary, hypo2Vec)

#### C1 analysis without hominins ####
#### Hypothesis 1: deciles 3-6; deciles 7-10 ####
testCR(lc_evolutionary_no_hominin, hypo1Vec)

#### Hypothesis 2: deciles 3-5; deciles 6-10 ####
testCR(lc_evolutionary_no_hominin, hypo2Vec)


# Rphylopars --------------------------------------------------------------

library(Rphylopars)
dat <- read.csv("Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv")
cn <- colnames(dat)
cn[1] <- "species"
colnames(dat) <- cn

tree <- ape::read.tree(file = "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt")
res <- phylopars(dat, tree)
saveRDS(res, paste0(input, "lc/phylopars/lc_dec3_10_phylopars.rds"))
phyloparsRes <-readRDS(paste0(input, "lc/phylopars/lc_dec3_10_phylopars.rds"))

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

lc_vcv_list <- readRDS(paste0(input, "lc/lc_dec3_10_vcv_extracted.RDS"))

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
lc_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_dec3_10_vcv_extracted.RDS")
lc_evolutionary <- lc_vcv_list$evolutionary
lc_vcv_list_no_hominins <- readRDS(paste0(input, "lc/lc_dec3_10_no_hominin_vcv_extracted.RDS"))
lc_evolutionary_no_hominin <- lc_vcv_list_no_hominins$evolutionary
ui2_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_dec3_10_no_pongo_vcv_extracted.RDS")
ui2_evolutionary <- ui2_vcv_list$evolutionary

library(evolqg)
library(parallel)


# helper functions --------------------------------------------------------
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


# run analysis ------------------------------------------------------------
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
