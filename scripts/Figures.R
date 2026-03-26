library(ape)
library(bayestestR)
library(dplyr)
library(MASS)
library(ggridges)
library(ggplot2)
library(ggtree)
library(gghalves)
library(overlapping)
library(parallel)
library(patchwork)
library(RColorBrewer)
library(tidyr)
library(tidyverse)

output <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/"

posterior <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8.tsv")
posterior <- posterior[round(0.1 * nrow(posterior)) : nrow(posterior), ] #apply burnin


# modern human line dat ---------------------------------------------------
lcDat <- read.csv("Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv")
mh <- filter(
  lcDat,
  genus == "Homo_sapiens"
)

df_long <- mh %>%
  mutate(id = row_number()) %>%
  pivot_longer(
    cols = c(Decile.3, Decile.4, Decile.5, Decile.6,
             Decile.7, Decile.8, Decile.9,
             `Buccal.decile.10..cervical.`),
    names_to  = "decile",
    values_to = "value"
  ) %>%
  mutate(
    decile_num = case_when(
      decile == "Decile.3"                    ~ "Decile 3",
      decile == "Decile.4"                    ~ "Decile 4",
      decile == "Decile.5"                    ~ "Decile 5",
      decile == "Decile.6"                    ~ "Decile 6",
      decile == "Decile.7"                    ~ "Decile 7",
      decile == "Decile.8"                    ~ "Decile 8",
      decile == "Decile.9"                    ~ "Decile 9",
      decile == "Buccal.decile.10..cervical." ~ "Decile 10"
    )
  )
df_long <- df_long %>%
  mutate(decile_num = factor(decile_num, levels = paste("Decile", 3:10)))

p1 <- ggplot(df_long, aes(x = decile_num, y = value, group = id)) +
  geom_line(alpha = 1.0, linewidth = 0.5, color = "black") +
  geom_point(alpha = 1.0, size = 0.8, color = "black") +
  labs(
    x = NULL,
    y = "Perikymata per millimeter",
  ) +
  theme_minimal(base_family = "Georgia") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p1
ggsave(paste0(output, "homoData.svg"), plot = p1, width = 7, height = 6)
  
# evo VCV, tree, tipmeans ---------------------------------------
map_estimate <- function(x) {
  as.numeric(bayestestR::map_estimate(x))
}

#evo VCV MAP
evo_vcv_cols <- paste0("evo_vcv_.", rep(0:7, each = 8), ".", rep(0:7, times = 8), ".")

evo_map <- dplyr::select(posterior, all_of(evo_vcv_cols)) |>
  summarise(across(everything(), map_estimate)) |>
  pivot_longer(everything(), names_to = "element", values_to = "map") |>
  mutate(
    row = as.integer(sub("evo_vcv_\\.(\\d+)\\.(\\d+)\\.", "\\1", element)) + 1,
    col = as.integer(sub("evo_vcv_\\.\\d+\\.(\\d+)\\.", "\\1", element)) + 1
  )

decile_labels <- paste0("Decile ", 3:10)

evo_map <- evo_map |>
  mutate(
    row_label = factor(decile_labels[row], levels = decile_labels),
    col_label = factor(decile_labels[col], levels = decile_labels)
  )

evoVCV <- ggplot(evo_map, aes(x = col_label, y = fct_rev(row_label), fill = map)) +
  geom_tile() +
  geom_text(aes(label = round(map, 2)), size = 3, color = "black") +
  scale_fill_gradient(
    low  = "white",
    high = "#a31e22"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "MAP"
  ) +
  theme_minimal(base_family = "Georgia") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
evoVCV

#tree
plottree <- ape::read.tree(file = "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt")
plottree$tip.label <- gsub("_", " ", plottree$tip.label)
plottree$tip.label <- gsub("Homo sapiens", "Modern humans", plottree$tip.label)
treeplot <- ggtree(plottree) + 
              geom_tiplab(aes(fontface = ifelse(label %in% c("Modern humans", "Neanderthal"), 2, 4)), family = "Georgia") +
              hexpand(0.55)
treeplot <- ggtree::rotate(treeplot, 12)
treeplot

#heatmap at tips
tip_order <- treeplot$data |>
  dplyr::filter(isTip) |>
  dplyr::arrange(y) |>
  dplyr::pull(label)

mean_map <- lapply(names(species_map), function(tip_label) {
  sp        <- unname(species_map[tip_label])
  mean_cols <- paste0(sp, "_mean_", 0:7)
  vals      <- dplyr::select(posterior, all_of(mean_cols)) |>
    summarise(across(everything(), map_estimate))
  data.frame(
    tip_label = tip_label,
    decile    = paste0("Decile ", 3:10),
    map       = as.numeric(vals)
  )
}) |> bind_rows() |>
  mutate(
    decile    = factor(decile, levels = paste0("Decile ", 3:10)),
    tip_label = factor(tip_label, levels = tip_order)
  )

italic_face <- ifelse(tip_order %in% c("Modern humans", "Neanderthal"), 
                      "plain", "italic")

tipplot <- ggplot(mean_map, aes(x = decile, y = tip_label, fill = map)) +
  geom_tile() +
  geom_text(aes(label = round(map, 1)),
            size   = 3.5,
            color  = "black",
            family = "Georgia") +
  scale_fill_gradient2(
    low  = "white",
    high = "#09539c",
    name = "MAP mean"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Georgia"),
    axis.text.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )
tipplot


p1 <- evoVCV + treeplot + tipplot
p1
ggsave(paste0(output, "treeMAP.svg"), plot = p1, width = 14, height = 6)

# Posterior predictive differences between modern humans and neand --------
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

n_samples <- nrow(posterior)
n_traits <- 8
trait_labels <- paste0("Decile ", 3:10)

# species_colors <- c(
#   "Modern humans"    = "#2166AC",
#   "Neanderthals"     = "#AB6621",
#   
#   "Pan paniscus"     = "#623A04",
#   "Pan troglodytes"  = "#042C62",
#   
#   "Gorilla beringei" = "#650109",
#   "Gorilla gorilla"  = "#01665E",
#   
#   "Pongo abelii"     = "#762A83",
#   "Pongo pygmaeus"   = "#37832A"
# )


# Generate posterior predictives
# hs_preds <- get_posterior_predictive(posterior, "Homo_sapiens", n_traits, n_samples)
# ne_preds <- get_posterior_predictive(posterior, "Neanderthal", n_traits, n_samples)
# pp_preds <- get_posterior_predictive(posterior, "Pan_paniscus", n_traits, n_samples)
# pt_preds <- get_posterior_predictive(posterior, "Pan_troglodytes", n_traits, n_samples)
# gb_preds <- get_posterior_predictive(posterior, "Gorilla_beringei", n_traits, n_samples)
# gg_preds <- get_posterior_predictive(posterior, "Gorilla_gorilla", n_traits, n_samples)
# pa_preds <- get_posterior_predictive(posterior, "Pongo_abelii", n_traits, n_samples)
# ppyg_preds <- get_posterior_predictive(posterior, "Pongo_pygmaeus", n_traits, n_samples)
# saveRDS(hs_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/hsPostPred.rds")
# saveRDS(ne_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/neanderthalPostPred.rds")
# saveRDS(pp_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/panpaniscusPostPred.rds")
# saveRDS(pt_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pantroglodytesPostPred.rds")
# saveRDS(gb_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/gorrillaberingeiPostPred.rds")
# saveRDS(gg_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/gorillagorillaPostPred.rds")
# saveRDS(pa_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pongoabeliiPostPred.rds")
# saveRDS(ppyg_preds, "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pongopygmaeusPostPred.rds")
hs_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/hsPostPred.rds")
ne_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/neanderthalPostPred.rds")
pp_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/panpaniscusPostPred.rds")
pt_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pantroglodytesPostPred.rds")
gb_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/gorrillaberingeiPostPred.rds")
gg_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/gorillagorillaPostPred.rds")
pa_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pongoabeliiPostPred.rds")
ppyg_preds <- read_rds("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/PosteriorPredictiveDraws/pongopygmaeusPostPred.rds")


plotRidgePlot <- function(pred1, pred2, specName1, specName2, plotName1, plotName2, color1, color2){
  recode_vec <- setNames(c(plotName1, plotName2), c(specName1, specName2))
  color_vec  <- setNames(c(color1, color2), c(plotName1, plotName2))
  
  plot_data <- bind_rows(pred1, pred2) |>
    pivot_longer(cols = all_of(trait_labels),
                 names_to = "trait",
                 values_to = "value") |>
    mutate(
      trait   = factor(trait, levels = trait_labels),
      species = recode(species, !!!recode_vec)
    )
  
  overlap_data <- plot_data |>
    group_by(trait) |>
    summarise(
      overlap = overlapping::overlap(
        list(
          value[species == plotName1],
          value[species == plotName2]
        )
      )$OV,
      .groups = "drop"
    ) |>
    mutate(
      trait = factor(trait, levels = trait_labels),
      label = paste0(round(overlap * 100, 1), "%")
    )
  
  pt <- ggplot(plot_data, aes(x = trait, y = value, fill = species)) +
    geom_half_violin(data = filter(plot_data, species == plotName1),
                     aes(fill = species),
                     alpha = 0.6, scale = "width", side = "l") +
    geom_half_violin(data = filter(plot_data, species == plotName2),
                     aes(fill = species),
                     alpha = 0.6, scale = "width", side = "r") +
    geom_half_boxplot(data = filter(plot_data, species == plotName1),
                     alpha = 0.6, scale = "width", side = "l", outlier.shape = NA) +
    geom_half_boxplot(data = filter(plot_data, species == plotName2),
                     alpha = 0.6, scale = "width", side = "r", outlier.shape = NA) +
    geom_text(data = overlap_data,
              aes(x = trait, y = Inf, label = label),
              inherit.aes = FALSE,
              vjust = 1.5, size = 3, color = "grey30") +
    scale_fill_manual(values = color_vec) +
    scale_y_continuous(limits = c(0, 50)) +
    labs(
      x    = "Decile",
      y    = "Perikymata count per millimeter",
      fill = "Species"
    ) +
    theme_minimal(base_family = "Georgia") +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
    )
  return(pt)
}

colors <- brewer.pal(8, "Paired")

species_colors <- c(
  "Modern humans"    = colors[1],
  "Neanderthals"     = colors[2],
  
  "Pan paniscus"     = colors[3],
  "Pan troglodytes"  = colors[4],
  
  "Gorilla beringei" = colors[5],
  "Gorilla gorilla"  = colors[6],
  
  "Pongo abelii"     = colors[7],
  "Pongo pygmaeus"   = colors[8]
)

homo <- plotRidgePlot(hs_preds, ne_preds, 
                      "Homo_sapiens", "Neanderthal", 
                      "Modern humans", "Neanderthals", 
                      species_colors["Modern humans"], species_colors["Neanderthals"])
homo

pan <- plotRidgePlot(pp_preds, pt_preds, 
                     "Pan_paniscus", "Pan_troglodytes", 
                     "Pan paniscus", "Pan troglodytes", 
                     species_colors["Pan paniscus"], species_colors["Pan troglodytes"])
pan

gorilla <- plotRidgePlot(gb_preds, gg_preds,
                         "Gorilla_beringei", "Gorilla_gorilla", 
                         "Gorilla beringei", "Gorilla gorilla", 
                         species_colors["Gorilla beringei"], species_colors["Gorilla gorilla"])
gorilla

pongo <- plotRidgePlot(ppyg_preds, pa_preds, 
                     "Pongo_abelii", "Pongo_pygmaeus", 
                     "Pongo abelii", "Pongo pygmaeus", 
                     species_colors["Pongo abelii"], species_colors["Pongo pygmaeus"])
pongo

combined <- homo + pan + gorilla + pongo
combined
ggsave(paste0(output, "postPred.svg"), plot = combined, width = 14, height = 14)



