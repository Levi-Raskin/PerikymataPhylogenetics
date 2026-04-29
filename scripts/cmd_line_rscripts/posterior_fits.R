args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

p <- 8

vcv_list <- readRDS(input)

posteriorFits <- list()

for (i in 1:length(vcv_list)) {
  post <- vcv_list[[i]]
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
  
  for (iter in 1:500) {
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

names(posteriorFits) <- names(vcv_list)
saveRDS(posteriorFits, output)

message("----------------------------------------")
message(paste0("Posterior IW fitting finished for: ", input))
message("----------------------------------------")