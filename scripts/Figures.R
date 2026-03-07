library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(MASS)
library(ggridges)
library(dplyr)
library(tidyr)
library(parallel)


plotdat <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8.tsv")

plotdat <- plotdat[round(0.1 * nrow(plotdat)) : nrow(plotdat), ] #apply burnin

# Posterior predictive differences between modern humans and neanderthals
n_samples <- nrow(plotdat)
n_traits <- 8
trait_labels <- paste0("Decile ", 3:10)

get_posterior_predictive <- function(posterior, species, n_traits, n_samples) {
  
  mean_cols <- paste0(species, "_mean_", 0:(n_traits - 1))
  mu_samples <- as.matrix(posterior[, mean_cols])
  
  vcv_cols <- outer(0:(n_traits - 1), 0:(n_traits - 1),
                    FUN = function(i, j) paste0(species, "_vcv_.", i, ".", j, "."))
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


# Generate posterior predictives for both species
hs_preds <- get_posterior_predictive(plotdat, "Homo_sapiens", n_traits, n_samples)
ne_preds <- get_posterior_predictive(plotdat, "Neanderthal", n_traits, n_samples)

# Combine and reshape to long format
plot_data <- bind_rows(hs_preds, ne_preds) |>
  pivot_longer(cols = all_of(trait_labels),
               names_to = "trait",
               values_to = "value") |>
  mutate(
    trait = factor(trait, levels = rev(trait_labels)),  # reverse so decile 1 is on top
    species = recode(species,
                     "Homo_sapiens" = "Modern Human",
                     "Neanderthal"  = "Neanderthal")
  )

# Plot
ggplot(plot_data, aes(x = value, y = trait, fill = species, color = species)) +
  geom_density_ridges(alpha = 0.4, scale = 0.9, rel_min_height = 0.01) +
  scale_fill_manual(values = c("Modern Human" = "#2166AC", "Neanderthal" = "#D6604D")) +
  scale_color_manual(values = c("Modern Human" = "#2166AC", "Neanderthal" = "#D6604D")) +
  labs(
    x = "Perikymata Count",
    y = NULL,
    fill = "Species",
    color = "Species"
  ) +
  theme_ridges(grid = FALSE) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10)
  )


# Estimated mean values ---------------------------------------------------

#to estimate MAP
posterior_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

plotMeans <- plotdat %>%
  select(matches("^intraspecificMean")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("taxon", "trait"),
    names_pattern = "intraspecificMean(.+)(\\d+)$"
  ) %>%
  mutate(trait = paste0("Decile ", as.integer(trait) + 3)) %>%
  group_by(taxon, trait) %>%
  summarise(mean_value = mean(value), .groups = "drop")

plotMeans$trait <- factor(plotMeans$trait, 
                          levels = c("Decile 3", 
                                     "Decile 4",
                                     "Decile 5",
                                     "Decile 6", 
                                     "Decile 7",
                                     "Decile 8", 
                                     "Decile 9",
                                     "Decile 10"
                                     ))

plotMeans$taxon <- gsub("_", plotMeans$taxon, replacement = " ")


plotmeansSister <- filter(plotMeans, taxon %in% c("Gorilla beringei", "Gorilla gorilla"))
p <- ggplot(plotmeansSister, aes(x = trait, y = mean_value, color = taxon, group = taxon)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Trait",
    y = "Mean inferred posterior mean value",
    color = "Taxon"
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(10, 17), breaks = 10:18) +
  theme_minimal()
p
ggsave("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/gorillasMeanDeciles.svg", p, width = 8, height = 5)

plotmeansSister <- filter(plotMeans, taxon %in% c("Pongo abelii", "Pongo pygmaeus"))
p <- ggplot(plotmeansSister, aes(x = trait, y = mean_value, color = taxon, group = taxon)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Trait",
    y = "Mean inferred posterior mean value",
    color = "Taxon"
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(10, 17), breaks = 10:18) +
  theme_minimal()
p
ggsave("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/pongosMeanDeciles.svg", p, width = 8, height = 5)


plotmeansSister <- filter(plotMeans, taxon %in% c("Pan troglodytes", "Pan paniscus"))
p <- ggplot(plotmeansSister, aes(x = trait, y = mean_value, color = taxon, group = taxon)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Trait",
    y = "Mean inferred posterior mean value",
    color = "Taxon"
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(10, 17), breaks = 10:18) +
  theme_minimal()
p
ggsave("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/pansMeanDeciles.svg", p, width = 8, height = 5)


plotmeansSister <- filter(plotMeans, taxon %in% c("Homo sapiens", "Neanderthal"))
p <- ggplot(plotmeansSister, aes(x = trait, y = mean_value, color = taxon, group = taxon)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Trait",
    y = "Mean inferred posterior mean value",
    color = "Taxon"
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(10, 17), breaks = 10:18) +
  theme_minimal()
p
ggsave("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/homininsMeanDeciles.svg", p, width = 8, height = 5)


p <- ggplot(plotMeans, aes(x = trait, y = mean_value, color = taxon, group = taxon)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Trait",
    y = "Mean inferred posterior mean value",
    color = "Taxon"
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(10, 17), breaks = 10:18) +
  theme_minimal()
p
ggsave("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/allMeanDeciles.svg", p, width = 8, height = 5)


# estimated variances -----------------------------------------------------

plotVCV <- plotdat %>%
  select(matches("^intraspecificVCV")) %>%
  select(matches("(\\d+)\\.\\1$")) %>%  # keep only diagonals (0.0, 1.1, 2.2, ...)
  pivot_longer(
    cols = everything(),
    names_to = c("taxon", "trait"),
    names_pattern = "intraspecificVCV(.+)(\\d+)\\.\\d+$"
  ) %>%
  mutate(trait = paste0("Decile ", as.integer(trait) + 3)) %>%
  group_by(taxon, trait) %>%
  summarise(mean_value = posterior_mode(value), .groups = "drop")

plotVCV$trait <- factor(plotVCV$trait, 
                        levels = c("Decile 3", 
                                   "Decile 4",
                                   "Decile 5",
                                   "Decile 6", 
                                   "Decile 7",
                                   "Decile 8", 
                                   "Decile 9",
                                   "Decile 10"))

plotVCV$taxon <- gsub("_", plotMeans$taxon, replacement = " ")

p <- ggplot(plotVCV, aes(x = trait, y = mean_value, color = taxon, group = taxon)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Trait",
    y = "Mean inferred variance",
    color = "Taxon"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal()
p
ggsave("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/allVariances.svg", p, width = 8, height = 5)

