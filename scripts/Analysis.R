library(coda)
library(dplyr)
library(MCMCpack)
library(overlapping)
library(parallel)

lc_posterior <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10.tsv"))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

lc_posterior_no_hominin <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10_no_hominin.tsv"))
lc_posterior_no_hominin  <- lc_posterior_no_hominin[round(0.1 * nrow(lc_posterior_no_hominin)) : nrow(lc_posterior_no_hominin), ] #apply burnin

ui2_posterior <- as.data.frame(fread("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/ui2/ui2_dec3_10_no_pongo.tsv"))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ] #apply burnin

# Functions ---------------------------------------------------------------
convertLatexTable <- function(vec){
  string <- ""
  for(i in vec){
    string <- paste0(string, round(i, 2), " & ")
  }
  print(string)
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

#GR
lc_gr <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10_ess_gelman_rubin.RDS")

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

ui2_gr <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/ui2/ui2_dec3_10_ess_gelman_rubin.RDS")

# KL divergence -------------------------------------------------
p <- 8
#### Lower Canine ####
lc_posteriorFits <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10_posterior_fits.RDS")

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
convertLatexTable(res)

for(i in 1:length(lc_posteriorFits)){
  print(paste(names(lc_posteriorFits)[i],round( lc_posteriorFits[[i]]$nu , 2)))
}

#### Lower Canine no hominins ####
lc_no_hominin_posteriorFits <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10_no_hominin_posterior_fits.RDS") 
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
convertLatexTable(res)

for(i in 1:length(lc_no_hominin_posteriorFits)){
  print(paste(names(lc_no_hominin_posteriorFits)[i],round( lc_no_hominin_posteriorFits[[i]]$nu , 2)))
}

#### Upper second incisor ####
ui2_posteriorFits <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/ui2/ui2_dec3_10_no_pongo_posterior_fits.RDS")
  
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
convertLatexTable(res)

for(i in 1:length(ui2_posteriorFits)){
  print(paste(names(ui2_posteriorFits)[i],round( ui2_posteriorFits[[i]]$nu , 2)))
}


#### Lower canine symmetrized KL divergence ####
lc_posteriorFits <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/lc_dec3_10_posterior_fits.RDS")

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

hs_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/hsPostPred.rds")
ne_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/neanderthalPostPred.rds")
pp_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/panpaniscusPostPred.rds")
pt_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/pantroglodytesPostPred.rds")
gb_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/gorrillaberingeiPostPred.rds")
gg_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/gorillagorillaPostPred.rds")
pa_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/pongoabeliiPostPred.rds")
ppyg_preds <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs/lc/posteriorPredictive/pongopygmaeusPostPred.rds")

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
lc_vcv_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_VCVs_extracted.rds")
lc_evolutionary <- lc_vcv_list$evolutionary
lc_vcv_list_no_hominins <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_no_hominin_VCVs_extracted.rds")
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
  
  cat("Posterior mean AVG ratio:", mean(avg_rat_results, na.rm = TRUE), "\n")
  cat("95% credible interval:", quantile(avg_rat_results, c(0.025, 0.975), na.rm = TRUE), "\n")
}

hypoMat1 <- matrix(0, 8, 8)
hypoMat1[1:4, 1:4] <- 1
hypoMat1[5:8, 5:8] <- 1

hypoMat2 <- matrix(0, 8, 8)
hypoMat2[1:4, 1:4] <- 1
hypoMat2[5:8, 5:8] <- 1

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

res <- phylopars(dat, tree)
