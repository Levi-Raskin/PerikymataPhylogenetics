library(ggplot2)
library(tidyverse)
library(RColorBrewer)

plotdat <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/LC_prelim/out0gelmanRubin2.tsv")
plotdat2 <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/LC_prelim/out0gelmanRubin0.tsv")
plotdat3 <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/LC_prelim/out0gelmanRubin1.tsv")
plotdat4 <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/LC_prelim/out0gelmanRubin3.tsv")

plotdat <- plotdat[round(0.1 * nrow(plotdat)) : nrow(plotdat), ] #apply burnin
plotdat2 <- plotdat[round(0.1 * nrow(plotdat)) : nrow(plotdat), ] #apply burnin
plotdat3 <- plotdat[round(0.1 * nrow(plotdat)) : nrow(plotdat), ] #apply burnin
plotdat4 <- plotdat[round(0.1 * nrow(plotdat)) : nrow(plotdat), ] #apply burnin

plotdat <- rbind(
  plotdat,
  plotdat2,
  plotdat3,
  plotdat4
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

