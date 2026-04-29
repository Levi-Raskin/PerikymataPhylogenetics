library(data.table)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

posterior <- as.data.frame(fread(input))
posterior <- posterior[round(0.1 * nrow(posterior)) : nrow(posterior), ] #apply burnin

extract_vcv <- function(row, prefix, p = 8) {
  mat <- matrix(NA, p, p)
  for (i in 0:(p-1)) {
    for (j in 0:(p-1)) {
      col <- paste0(prefix, "vcv_(", i, ",", j, ")")
      if (col %in% names(row)) {
        mat[i+1, j+1] <- row[[col]]
      }
    }
  }
  mat
}

species_list <- c(
  "evolutionary",
  "Pongo_abelii",
  "Pongo_pygmaeus",
  "Pan_troglodytes",
  "Pan_paniscus",
  "Gorilla_beringei",
  "Gorilla_gorilla",
  "Homo_sapiens",
  "Neanderthal"
)

prefix_map <- c(
  evolutionary    = "evo_",
  Pongo_abelii    = "Pongo_abelii_",
  Pongo_pygmaeus  = "Pongo_pygmaeus_",
  Pan_troglodytes = "Pan_troglodytes_",
  Pan_paniscus    = "Pan_paniscus_",
  Gorilla_beringei = "Gorilla_beringei_",
  Gorilla_gorilla  = "Gorilla_gorilla_",
  Homo_sapiens     = "Homo_sapiens_",
  Neanderthal      = "Neanderthal_"
)

present_species <- species_list[sapply(species_list, function(sp) {
  any(startsWith(colnames(posterior), prefix_map[sp]))
})]

message("Species found in posterior: ", paste(present_species, collapse = ", "))

n_samples <- nrow(posterior)

vcv_list <- lapply(present_species, function(sp) {
  prefix <- prefix_map[sp]
  mclapply(seq_len(n_samples), function(i) {
    extract_vcv(posterior[i, ], prefix, p = 8)
  },
  mc.cores = detectCores() - 1)
})

names(vcv_list) <- present_species
saveRDS(vcv_list, output)

message("----------------------------------------")
message(paste0("VCV list processing finished for: ", input))
message("----------------------------------------")