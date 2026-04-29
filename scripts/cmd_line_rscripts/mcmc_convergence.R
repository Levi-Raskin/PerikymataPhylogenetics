library(data.table)
library(parallel)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

# read posteriors and form into list
post_list <- list()
for(i in 1:4){
  post <- as.data.frame(fread(paste0(input, i, ".tsv")))
  post <- post[round(0.1 * nrow(post)) : nrow(post), 3:ncol(post)] #apply burnin, remove cycle and lnl 
  post_list[[i]] <- post
}

#calculate rank normalized r hat for each parameter
nparm <- ncol(post_list[[1]])
ncyc <- nrow(post_list[[1]])
nchains = 4
parmnames <- names(post_list[[1]])

calcRankNormRhat <- function(idx){
  mat <- vapply(post_list, function(chain) chain[[idx]], numeric(ncyc))
  posterior::rhat(mat)
}

res <- mclapply(1:nparm,
                calcRankNormRhat,
                mc.cores = detectCores() - 1)
names(res) = parmnames

message("----------------------------------------")
message(paste0("Rank normalized r hat finished for: ", input))
message(paste0("Mean r hat: ", mean(unlist(res))))
message(paste0("Max r hat: ", max(unlist(res))))
message("----------------------------------------")
saveRDS(res, output)
