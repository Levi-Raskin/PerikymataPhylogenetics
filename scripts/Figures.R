library(ape)
library(bayestestR)
library(data.table)
library(dplyr)
library(latex2exp)
library(MASS)
library(ggridges)
library(ggplot2)
library(ggtree)
library(gghalves)
library(overlapping)
library(parallel)
library(patchwork)
library(RColorBrewer)
library(RiemBase)
library(tidyr)
library(tidyverse)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::mutate)

input <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_MCMC/"
output <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/"

lc_posterior <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10.tsv")))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

# lc_posterior_no_hominin <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_no_hominin.tsv")))
# lc_posterior_no_hominin  <- lc_posterior_no_hominin[round(0.1 * nrow(lc_posterior_no_hominin)) : nrow(lc_posterior_no_hominin), ] #apply burnin

lc_posterior_species_means <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_species_means.tsv")))
lc_posterior_species_means  <- lc_posterior_species_means[round(0.1 * nrow(lc_posterior_species_means)) : nrow(lc_posterior_species_means), ] #apply burnin

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ] #apply burnin

ui2_posterior_species_means <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_species_means.tsv")))
ui2_posterior_species_means  <- ui2_posterior_species_means[round(0.1 * nrow(ui2_posterior_species_means)) : nrow(ui2_posterior_species_means), ] #apply burnin


# modern human line drawing ---------------------------------------------------
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

lc_vcv_list                <- readRDS(paste0(input, "lc/lc_dec3_10_vcv_extracted.RDS"))
ui2_vcv_list                <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_vcv_extracted.RDS"))
lc_vcv_list_species_means   <- readRDS(paste0(input, "lc/lc_dec3_10_vcv_extracted_species_means.RDS"))
ui2_vcv_list_species_means  <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_vcv_extracted_species_means.RDS"))

vcv_mean_estimate <- function(vcv_list) {
  data_list <- riemfactory(vcv_list, name = "spd")
  x <- rbase.mean(data_list)
  return(x$x)
}

decile_labels <- paste0("Decile ", 3:10)

matrix_to_long <- function(mat) {
  df <- as.data.frame(mat)
  colnames(df) <- 0:(ncol(mat) - 1)
  df$row <- 0:(nrow(mat) - 1)
  
  tidyr::pivot_longer(df, -row, names_to = "col", values_to = "mean_val") |>
    mutate(
      col       = as.integer(col),
      row_label = factor(decile_labels[row + 1], levels = decile_labels),
      col_label = factor(decile_labels[col + 1], levels = decile_labels),
      is_diag   = row_label == col_label
    )
}

# Compute each group's Fréchet-mean evolutionary VCV matrix once
lc_mean     <- vcv_mean_estimate(lc_vcv_list$evolutionary)
ui2_mean    <- vcv_mean_estimate(ui2_vcv_list$evolutionary)
lc_mean_sp  <- vcv_mean_estimate(lc_vcv_list_species_means$evolutionary)
ui2_mean_sp <- vcv_mean_estimate(ui2_vcv_list_species_means$evolutionary)

shared_max <- max(lc_mean, ui2_mean, lc_mean_sp, ui2_mean_sp)
shared_min <- min(lc_mean, ui2_mean, lc_mean_sp, ui2_mean_sp)

# Helper to build a heatmap so we're not retyping the whole ggplot call four times
build_evoVCV <- function(mean_mat) {
  evo_map <- matrix_to_long(mean_mat)
  
  ggplot(evo_map, aes(x = col_label, y = fct_rev(row_label), fill = mean_val)) +
    geom_tile() +
    geom_tile(data = subset(evo_map, is_diag), color = "black", linewidth = 1.5, fill = NA) +
    geom_text(aes(label = round(mean_val, 2)), size = 3, color = "black") +
    scale_fill_gradient(
      low    = "white",
      high   = "#a31e22",
      limits = c(shared_min, shared_max)
    ) +
    labs(x = NULL, y = NULL, fill = "Mean") +
    theme_minimal(base_family = "Georgia") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

### LC full ###
evoVCV_lc <- build_evoVCV(lc_mean)
evoVCV_lc

#tree
plottree_lc <- ape::read.tree(file = "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt")
plottree_lc$tip.label <- gsub("_", " ", plottree_lc$tip.label)
plottree_lc$tip.label <- gsub("Homo sapiens", "Modern humans", plottree_lc$tip.label)
treeplot_lc <- ggtree(plottree_lc) + 
  geom_tiplab(aes(fontface = ifelse(label %in% c("Modern humans", "Neanderthal"), 2, 4)), family = "Georgia") +
  hexpand(0.55)
treeplot_lc <- ggtree::rotate(treeplot_lc, 12)
treeplot_lc <- ggtree::rotate(treeplot_lc, 14)
treeplot_lc

lcHominins <- (evoVCV_lc + treeplot_lc) + plot_layout(widths = c(1, 1))
lcHominins

### LC wo intra (species means) ###
evoVCV_lc_sp <- build_evoVCV(lc_mean_sp)
evoVCV_lc_sp

treeplot_lc_sp <- treeplot_lc  # same tree/rotation as above

lcHominins_species <- (evoVCV_lc_sp + treeplot_lc_sp) + plot_layout(widths = c(1, 1))
lcHominins_species

### UI2 ###
evoVCV_ui2 <- build_evoVCV(ui2_mean)
evoVCV_ui2

plottree_ui2 <- ape::read.tree(file = "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt")
plottree_ui2 <- keep.tip(plottree_ui2, c("Neanderthal", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"))
plottree_ui2$tip.label <- gsub("_", " ", plottree_ui2$tip.label)
treeplot_ui2 <- ggtree(plottree_ui2) + 
  geom_tiplab(aes(fontface = ifelse(label %in% c("Modern humans", "Neanderthal"), 2, 4)), family = "Georgia") +
  hexpand(0.55)
treeplot_ui2 <- ggtree::rotate(treeplot_ui2, 5)
treeplot_ui2

ui2 <- (evoVCV_ui2 + treeplot_ui2) + plot_layout(widths = c(1, 1))
ui2

### UI2 species means ###
evoVCV_ui2_sp <- build_evoVCV(ui2_mean_sp)
evoVCV_ui2_sp

treeplot_ui2_sp <- treeplot_ui2  # same pruned tree/rotation as above

ui2_species_means <- (evoVCV_ui2_sp + treeplot_ui2_sp) + plot_layout(widths = c(1, 1))
ui2_species_means

### Combine as one flat grid — avoids nested-patchwork width bugs ###
combined <- wrap_plots(
  evoVCV_lc,     treeplot_lc,
  evoVCV_lc_sp,  treeplot_lc_sp,
  evoVCV_ui2,    treeplot_ui2,
  evoVCV_ui2_sp, treeplot_ui2_sp,
  ncol = 4, nrow = 2,
  widths  = rep(1, 4),
  heights = c(1, 1)
)
combined

ggsave(paste0(output, "treeMAPCombined.svg"), plot = combined, width = 28, height = 12)


# SOM mean intraspecific VCV matrices -------------------------------------
# SOM mean intraspecific VCV matrices -------------------------------------
lc_vcv_list  <- readRDS(paste0(input, "lc/lc_dec3_10_vcv_extracted.RDS"))
ui2_vcv_list <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_vcv_extracted.RDS"))

vcv_mean_estimate <- function(vcv_list) {
  data_list <- riemfactory(vcv_list, name = "spd")
  x <- rbase.mean(data_list)
  return(x$x)
}

decile_labels <- paste0("Decile ", 3:10)

matrix_to_long <- function(mat) {
  df <- as.data.frame(mat)
  colnames(df) <- 0:(ncol(mat) - 1)
  df$row <- 0:(nrow(mat) - 1)
  
  tidyr::pivot_longer(df, -row, names_to = "col", values_to = "mean_val") |>
    mutate(
      col       = as.integer(col),
      row_label = factor(decile_labels[row + 1], levels = decile_labels),
      col_label = factor(decile_labels[col + 1], levels = decile_labels),
      is_diag   = row_label == col_label
    )
}

build_evoVCV <- function(mean_mat) {
  evo_map <- matrix_to_long(mean_mat)
  
  ggplot(evo_map, aes(x = col_label, y = fct_rev(row_label), fill = mean_val)) +
    geom_tile() +
    geom_tile(data = subset(evo_map, is_diag), color = "black", linewidth = 1.5, fill = NA) +
    geom_text(aes(label = round(mean_val, 2)), size = 3, color = "black") +
    scale_fill_gradient(
      low    = "white",
      high   = "#a31e22",
      limits = c(shared_min, shared_max)
    ) +
    labs(x = NULL, y = NULL, fill = "Mean") +
    theme_minimal(base_family = "Georgia") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

intraTaxonOrder <- c(
  "Homo_sapiens",
  "Neanderthal",
  "Pan_paniscus",
  "Pan_troglodytes",
  "Gorilla_beringei",
  "Gorilla_gorilla",
  "Pongo_abelii",
  "Pongo_pygmaeus"
)

### LC ###
lc_vcv_means <- list()
for (nm in intraTaxonOrder) {
  if (nm %in% names(lc_vcv_list)) {
    lc_vcv_means[[nm]] <- vcv_mean_estimate(lc_vcv_list[[nm]])
  }
}
lc_shared_max <- max(unlist(lc_vcv_means))
lc_shared_min <- min(unlist(lc_vcv_means))
shared_max <- lc_shared_max
shared_min <- lc_shared_min
lc_vcv_plots <- Map(function(mat, nm) build_evoVCV(mat) + ggtitle(gsub("_", " ", nm)) + theme(plot.title = element_text(hjust = 0.5)), lc_vcv_means, names(lc_vcv_means))

### UI2 ###
ui2_vcv_means <- list()
for (nm in intraTaxonOrder) {
  if (nm %in% names(ui2_vcv_list)) {
    ui2_vcv_means[[nm]] <- vcv_mean_estimate(ui2_vcv_list[[nm]])
  }
}
ui2_shared_max <- max(unlist(ui2_vcv_means))
ui2_shared_min <- min(unlist(ui2_vcv_means))
shared_max <- ui2_shared_max
shared_min <- ui2_shared_min
ui2_vcv_plots <- Map(function(mat, nm) build_evoVCV(mat) + ggtitle(gsub("_", " ", nm)) + theme(plot.title = element_text(hjust = 0.5)), ui2_vcv_means, names(ui2_vcv_means))

### Combine ###
n_col <- max(length(lc_vcv_plots), length(ui2_vcv_plots))
lc_vcv_plots  <- c(lc_vcv_plots,  replicate(n_col - length(lc_vcv_plots),  plot_spacer(), simplify = FALSE))
ui2_vcv_plots <- c(ui2_vcv_plots, replicate(n_col - length(ui2_vcv_plots), plot_spacer(), simplify = FALSE))

intraVCV_combined <- wrap_plots(c(lc_vcv_plots, ui2_vcv_plots), ncol = n_col, nrow = 2)
intraVCV_combined

ggsave(paste0(output, "intraspecificVCV_SOM.svg"), plot = intraVCV_combined, width = 5 * n_col, height = 10)

# UI2 LC combined ---------------------------------------------------------
hs_preds_lc   <- read_rds(paste0(input, "/lc/posteriorPredictive/hsPostPred.rds"))
ne_preds_lc   <- read_rds(paste0(input, "lc/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds_lc   <- read_rds(paste0(input, "lc/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds_lc   <- read_rds(paste0(input, "lc/posteriorPredictive/pantroglodytesPostPred.rds"))
gb_preds_lc   <- read_rds(paste0(input, "lc/posteriorPredictive/gorrillaberingeiPostPred.rds"))
gg_preds_lc   <- read_rds(paste0(input, "lc/posteriorPredictive/gorillagorillaPostPred.rds"))
pa_preds_lc   <- read_rds(paste0(input, "lc/posteriorPredictive/pongoabeliiPostPred.rds"))
ppyg_preds_lc <- read_rds(paste0(input, "lc/posteriorPredictive/pongopygmaeusPostPred.rds"))

hs_preds_ui2 <- read_rds(paste0(input, "ui2/posteriorPredictive/hsPostPred.rds"))
ne_preds_ui2 <- read_rds(paste0(input, "ui2/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds_ui2 <- read_rds(paste0(input, "ui2/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds_ui2 <- read_rds(paste0(input, "ui2/posteriorPredictive/pantroglodytesPostPred.rds"))

trait_labels <- paste0("Decile ", 3:10)

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

plotRidgePlotCombined <- function(pred1_lc, pred2_lc, pred1_ui2, pred2_ui2,
                                  specName1, specName2,
                                  plotName1, plotName2,
                                  color1, color2) {
  
  recode_vec <- setNames(c(plotName1, plotName2), c(specName1, specName2))
  color_vec  <- setNames(c(color1, color2), c(plotName1, plotName2))
  
  prep_data <- function(pred1, pred2, tooth_type) {
    bind_rows(pred1, pred2) |>
      pivot_longer(cols = all_of(trait_labels),
                   names_to = "trait",
                   values_to = "value") |>
      mutate(
        trait      = factor(trait, levels = trait_labels),
        species    = recode(species, !!!recode_vec),
        tooth_type = tooth_type
      )
  }
  
  ordered_levels <- c(
    paste(trait_labels, "C",  sep = " — "),
    paste(trait_labels, "I2", sep = " — ")
  )
  
  # Build display labels: "Decile 3 — $\bar{C}$", "Decile 3 — $I^2$", etc.
  x_labels <- setNames(
    c(sapply(trait_labels, function(d) TeX(paste0(d, " — $C_1$"))),
      sapply(trait_labels, function(d) TeX(paste0(d, " — $I^2$")))),
    ordered_levels
  )
  
  plot_data <- bind_rows(
    prep_data(pred1_lc,  pred2_lc,  "C"),
    prep_data(pred1_ui2, pred2_ui2, "I2")
  ) |>
    mutate(
      tooth_type  = factor(tooth_type, levels = c("C", "I2")),
      trait_tooth = factor(paste(trait, tooth_type, sep = " — "), levels = ordered_levels)
    )
  
  n_lc      <- length(trait_labels)
  band_data <- tibble(xmin = n_lc + 0.5, xmax = n_lc * 2 + 0.5)
  
  overlap_data <- plot_data |>
    group_by(trait_tooth) |>
    group_modify(~{
      df <- .x
      ov <- overlapping::overlap(
        list(
          df$value[df$species == plotName1],
          df$value[df$species == plotName2]
        )
      )$OV
      tibble(overlap = ov)
    }) |>
    ungroup() |>
    mutate(
      trait_tooth = factor(trait_tooth, levels = ordered_levels),
      label       = paste0(round(overlap * 100, 1), "%")
    )
  
  fill_vec <- setNames(
    c(color1, color2),
    c(plotName1, plotName2)
  )
  
  pt <- ggplot(plot_data, aes(x = trait_tooth, y = value, fill = species)) +
    geom_rect(
      data        = band_data,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill        = "grey90", alpha = 0.5
    ) +
    geom_half_violin(
      data  = filter(plot_data, species == plotName1),
      aes(fill = species),
      alpha = 0.7, scale = "width", side = "l"
    ) +
    geom_half_violin(
      data  = filter(plot_data, species == plotName2),
      aes(fill = species),
      alpha = 0.7, scale = "width", side = "r"
    ) +
    geom_half_boxplot(
      data  = filter(plot_data, species == plotName1),
      aes(fill = species),
      alpha = 0.7, scale = "width", side = "l", outlier.shape = NA
    ) +
    geom_half_boxplot(
      data  = filter(plot_data, species == plotName2),
      aes(fill = species),
      alpha = 0.7, scale = "width", side = "r", outlier.shape = NA
    ) +
    geom_text(
      data        = overlap_data,
      aes(x = trait_tooth, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.5, size = 5, color = "grey30"
    ) +
    scale_fill_manual(values = fill_vec) +
    scale_y_continuous(limits = c(0, 40)) +
    scale_x_discrete(labels = x_labels) +
    labs(
      x    = "Decile — Tooth type",
      y    = "Perikymata count per millimeter",
      fill = "Species"
    ) +
    theme_minimal(base_family = "Georgia") +
    theme(
      legend.position = "right",
      axis.text.x     = element_text(size = 12, angle = 45, hjust = 1)
    )
  
  return(pt)
}

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
    group_modify(~{
      df <- .x
      
      ov <- overlapping::overlap(
        list(
          df$value[df$species == plotName1],
          df$value[df$species == plotName2]
        )
      )$OV
      
      tibble(overlap = ov)
    }) |>
    ungroup() |>
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
              vjust = 1.5, size = 5, color = "grey30") +
    scale_fill_manual(values = color_vec) +
    scale_y_continuous(limits = c(0, 40)) +
    labs(
      x    = "Decile",
      y    = "Perikymata count per millimeter",
      fill = "Species"
    ) +
    theme_minimal(base_family = "Georgia") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
    )
  return(pt)
}
homo <- plotRidgePlotCombined(
  hs_preds_lc, ne_preds_lc, hs_preds_ui2, ne_preds_ui2,
  "Homo_sapiens", "Neanderthal",
  "Modern humans", "Neanderthals",
  species_colors["Modern humans"], species_colors["Neanderthals"]
)
homo

pan <- plotRidgePlotCombined(
  pp_preds_lc, pt_preds_lc, pp_preds_ui2, pt_preds_ui2,
  "Pan_paniscus", "Pan_troglodytes",
  "Pan paniscus", "Pan troglodytes",
  species_colors["Pan paniscus"], species_colors["Pan troglodytes"]
)
pan

# LC-only plots for gorilla and pongo (no UI2 data)
gorilla <- plotRidgePlot(
  gb_preds_lc, gg_preds_lc,
  "Gorilla_beringei", "Gorilla_gorilla",
  "Gorilla beringei", "Gorilla gorilla",
  species_colors["Gorilla beringei"], species_colors["Gorilla gorilla"]
)

pongo <- plotRidgePlot(
  ppyg_preds_lc, pa_preds_lc,
  "Pongo_abelii", "Pongo_pygmaeus",
  "Pongo abelii", "Pongo pygmaeus",
  species_colors["Pongo abelii"], species_colors["Pongo pygmaeus"]
)

combined <- homo /  pan / (gorilla + pongo)
combined

ggsave(paste0(output, "postPredCombined.svg"), plot = combined, width = 14, height = 20)

# Table 1: Posterior predictive means and variances -----------------------

writeTableOne <- function(pred){
  string <- ""
  for(i in 1:8){
    mean <- mean(pred[,i])
    var <- var(pred[,i])
    string <- paste0(string, round(mean, 2), " (" , round(var, 2), ")")
    if(i != 8){
      string <- paste0(string, " & ") 
    }
  }
  print(string)
}

#### LC
hs_preds <- read_rds(paste0(input, "lc/posteriorPredictive/hsPostPred.rds"))
ne_preds <- read_rds(paste0(input, "lc/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds <- read_rds(paste0(input, "lc/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds <- read_rds(paste0(input, "lc/posteriorPredictive/pantroglodytesPostPred.rds"))
gb_preds <- read_rds(paste0(input, "lc/posteriorPredictive/gorrillaberingeiPostPred.rds"))
gg_preds <- read_rds(paste0(input, "lc/posteriorPredictive/gorillagorillaPostPred.rds"))
pa_preds <- read_rds(paste0(input, "lc/posteriorPredictive/pongoabeliiPostPred.rds"))
ppyg_preds <- read_rds(paste0(input, "lc/posteriorPredictive/pongopygmaeusPostPred.rds"))

writeTableOne(hs_preds)
writeTableOne(ne_preds)
writeTableOne(pp_preds)
writeTableOne(pt_preds)
writeTableOne(gb_preds)
writeTableOne(gg_preds)
writeTableOne(pa_preds)
writeTableOne(ppyg_preds)

#### ui2
hs_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/hsPostPred.rds"))
ne_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/pantroglodytesPostPred.rds"))
writeTableOne(hs_preds)
writeTableOne(ne_preds)
writeTableOne(pp_preds)
writeTableOne(pt_preds)

# intraspecific means vs. MLE ------------------------------------------

plot_species_posteriors <- function(lc_posterior, lc_mle, species, bins = 500) {
  lc_mle <- lc_mle %>%
    dplyr::rename(any_of(c("Decile.10" = "Buccal.decile.10..cervical.")))
  decile_map <- c(
    "mean_0" = "Decile.3",
    "mean_1" = "Decile.4",
    "mean_2" = "Decile.5",
    "mean_3" = "Decile.6",
    "mean_4" = "Decile.7",
    "mean_5" = "Decile.8",
    "mean_6" = "Decile.9",
    "mean_7" = "Decile.10"
  )
  
  decile_levels <- paste0("Decile.", 3:10)
  cols          <- paste0(species, "_mean_", 0:7)
  
  missing_cols <- setdiff(cols, colnames(lc_posterior))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Species '%s' not found in lc_posterior. Missing columns: %s",
      species, paste(missing_cols, collapse = ", ")
    ))
  }
  
  if (all(!(species %in% lc_mle$genus))) {
    stop(sprintf(
      "Species '%s' not found in lc_mle$genus. Available species: %s",
      species, paste(unique(lc_mle$genus), collapse = ", ")
    ))
  }
  
  sp_long <- lc_posterior[, cols, drop = FALSE] %>%
    setNames(names(decile_map)) %>%
    pivot_longer(cols     = everything(),
                 names_to  = "mean_col",
                 values_to = "value") %>%
    mutate(decile = factor(unname(decile_map[mean_col]), levels = decile_levels))  # fixed
  
  sp_quantiles <- data.frame(
    decile = factor(decile_levels, levels = decile_levels),
    q_lo   = tapply(sp_long$value, sp_long$decile, quantile, 0.025),
    q_hi   = tapply(sp_long$value, sp_long$decile, quantile, 0.975)
  )
  
  sp_long <- sp_long %>%
    left_join(sp_quantiles, by = "decile") %>%
    mutate(
      fill_color = case_when(
        value <= q_lo ~ "tail",
        value >= q_hi ~ "tail",
        TRUE          ~ "middle"
      )
    )
  
  mle_sp <- lc_mle %>%
    filter(genus == species) %>%
    pivot_longer(cols      = starts_with("Decile."),
                 names_to  = "decile",
                 values_to = "mle_mean") %>%
    mutate(decile = factor(decile, levels = decile_levels))
  
  species_label <- gsub("_", " ", species)
  
  ggplot(sp_long, aes(x = value, fill = fill_color)) +
    geom_histogram(bins = bins, color = NA) +
    geom_vline(
      data      = mle_sp,
      aes(xintercept = mle_mean),
      color     = "blue",
      linewidth = 0.8,
      linetype  = "solid"
    ) +
    scale_fill_manual(
      values = c("tail" = "red", "middle" = "black"),
      labels = c("tail" = "Lower/Upper 2.5%", "middle" = "Middle 95%"),
      name   = NULL
    ) +
    facet_wrap(
      ~ decile, ncol = 4, scales = "free",
      labeller = labeller(decile = setNames(paste0("Decile ", 3:10), decile_levels))
    ) +
    labs(
      x       = "Inferred mean perikymata per millimeter",
      y       = "Posterior sample count",
      title   = bquote(italic(.(species_label)))
    ) +
    theme_minimal(base_family = "Georgia") +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}

species<- c(
  "Homo_sapiens",
  "Neanderthal",
  "Pan_paniscus",
  "Pan_troglodytes",
  "Gorilla_beringei",
  "Gorilla_gorilla",
  "Pongo_abelii",
  "Pongo_pygmaeus"
)

lc_mle <- read.csv("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10_species_means.csv")

for(i in species){
  p <- plot_species_posteriors(lc_posterior, lc_mle, i)
  print(p)
  ggsave(
    paste0(output, "/meanPosteriorHists/", i, ".pdf"), 
    plot = p, 
    width = 10, height = 8,
    device = cairo_pdf
  )
}

species<- c(
  "Homo_sapiens",
  "Neanderthal",
  "Pan_paniscus",
  "Pan_troglodytes"
)

ui2_mle <- read.csv("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo_species_means.csv")

for(i in species){
  p <- plot_species_posteriors(ui2_posterior, ui2_mle, i)
  print(p)
  ggsave(
    paste0(output, "/meanPosteriorHists/", i, "_ui2.pdf"), 
    plot = p, 
    width = 10, height = 8,
    device = cairo_pdf
  )
}


# Phylopars AIRM histograms -----------------------------------------------

lc_vcv_list <- readRDS(paste0(input, "lc/lc_dec3_10_vcv_extracted.RDS"))

all_AIRM_dat <- list()
for (i in 1:length(lc_vcv_list)) {
  name <- names(lc_vcv_list)[i]
  AIRM_dat <- readRDS(paste0(
    paste0(input, "lc/phylopars/",
           name,
           "_AIRM_distances.rds"
    )))
  all_AIRM_dat[[i]] <- data.frame(
    value = as.numeric(AIRM_dat),
    group = name
  )
}

colors <- brewer.pal(9, "Spectral")
combined_dat <- do.call(rbind, all_AIRM_dat)

p1 <- ggplot(combined_dat, aes(x = value, fill = group)) +
  geom_histogram(bins = 500, color = NaN, alpha = 0.7,
                 position = "identity") +
  scale_fill_manual(values = colors) +
  labs(
    x    = "AIRM distance",
    y    = "Count",
    fill = "VCV type"
  ) +
  theme_minimal(base_family = "Georgia") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )
p1
ggsave(
  paste0(output, "AIRM_distances_lc.pdf"), 
  plot = p1, 
  width = 10, height = 6,
  device = cairo_pdf
)


# posterior predictive check ----------------------------------------------
lc_data <- read.csv("Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv")
ui2_data <- read.csv("Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo.csv")

plot_ppc <- function(ppc_means, ppc_sds, observed_data, species, bins = 200) {
  library(patchwork)
  
  observed_data <- observed_data %>%
    dplyr::rename(any_of(c("Decile.10" = "Buccal.decile.10..cervical.")))
  decile_levels <- paste0("Decile.", 3:10)
  decile_labels <- setNames(paste0("Decile ", 3:10), decile_levels)
  
  obs_sp <- observed_data %>%
    filter(genus == species) %>%
    summarise(across(starts_with("Decile."),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd   = ~sd(.x,   na.rm = TRUE)))) %>%
    pivot_longer(everything(),
                 names_to  = c("decile", ".value"),
                 names_pattern = "(Decile\\.\\d+)_(mean|sd)") %>%
    mutate(decile = factor(decile, levels = decile_levels))
  
  ppc_means_long <- as.data.frame(ppc_means) %>%
    setNames(decile_levels) %>%
    pivot_longer(everything(), names_to = "decile", values_to = "value") %>%
    mutate(decile = factor(decile, levels = decile_levels))
  
  ppc_sds_long <- as.data.frame(ppc_sds) %>%
    setNames(decile_levels) %>%
    pivot_longer(everything(), names_to = "decile", values_to = "value") %>%
    mutate(decile = factor(decile, levels = decile_levels))
  
  obs_means <- obs_sp %>% transmute(decile, obs_value = mean)
  obs_sds   <- obs_sp %>% transmute(decile, obs_value = sd)
  
  pval_means <- ppc_means_long %>%
    left_join(obs_means, by = "decile") %>%
    group_by(decile) %>%
    summarise(pval = mean(value >= obs_value), .groups = "drop")
  
  pval_sds <- ppc_sds_long %>%
    left_join(obs_sds, by = "decile") %>%
    group_by(decile) %>%
    summarise(pval = mean(value >= obs_value), .groups = "drop")
  
  mean_labels <- setNames(
    paste0(decile_labels, "\np = ", round(pval_means$pval, 2)),
    decile_levels
  )
  sd_labels <- setNames(
    paste0(decile_labels, "\np = ", round(pval_sds$pval, 2)),
    decile_levels
  )
  
  species_label <- gsub("_", " ", species)
  
  p_mean <- ggplot(ppc_means_long, aes(x = value)) +
    geom_histogram(bins = bins, color = NA, fill = "black") +
    geom_vline(data = obs_means, aes(xintercept = obs_value),
               color = "blue", linewidth = 0.8) +
    facet_wrap(~ decile, nrow = 1, scales = "free",
               labeller = labeller(decile = mean_labels)) +
    scale_x_continuous(
      breaks = function(x) {
        brks <- scales::pretty_breaks(n = 4)(x)
        brks[brks %% 1 == 0]
      }
    ) +
    labs(x = "Posterior predictive mean", y = "Count",
         title = bquote(italic(.(species_label)))) +
    theme_minimal(base_family = "Georgia") +
    theme(panel.grid.minor = element_blank())
  
  p_sd <- ggplot(ppc_sds_long, aes(x = value)) +
    geom_histogram(bins = bins, color = NA, fill = "black") +
    geom_vline(data = obs_sds, aes(xintercept = obs_value),
               color = "blue", linewidth = 0.8) +
    facet_wrap(~ decile, nrow = 1, scales = "free",
               labeller = labeller(decile = sd_labels)) +
    scale_x_continuous(
      breaks = function(x) {
        brks <- scales::pretty_breaks(n = 6)(x)
        brks[brks %% 1 == 0]
      }
    ) +
    labs(x = "Posterior predictive SD", y = "Count") +
    theme_minimal(base_family = "Georgia") +
    theme(panel.grid.minor = element_blank())
  
  p_mean / p_sd
}

calc_ppc_pvalues <- function(ppc_means, ppc_sds, observed_data, species) {
  
  observed_data <- observed_data %>%
    rename(any_of(c("Decile.10" = "Buccal.decile.10..cervical.")))
  decile_levels <- paste0("Decile.", 3:10)
  
  obs_sp <- observed_data %>%
    filter(genus == species) %>%
    summarise(across(starts_with("Decile."),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd   = ~sd(.x,   na.rm = TRUE)))) %>%
    pivot_longer(everything(),
                 names_to  = c("decile", ".value"),
                 names_pattern = "(Decile\\.\\d+)_(mean|sd)") %>%
    mutate(decile = factor(decile, levels = decile_levels))
  
  obs_means <- obs_sp %>% transmute(decile, obs_value = mean)
  obs_sds   <- obs_sp %>% transmute(decile, obs_value = sd)
  
  ppc_means_long <- as.data.frame(ppc_means) %>%
    setNames(decile_levels) %>%
    pivot_longer(everything(), names_to = "decile", values_to = "value") %>%
    mutate(decile = factor(decile, levels = decile_levels))
  
  ppc_sds_long <- as.data.frame(ppc_sds) %>%
    setNames(decile_levels) %>%
    pivot_longer(everything(), names_to = "decile", values_to = "value") %>%
    mutate(decile = factor(decile, levels = decile_levels))
  
  pval_means <- ppc_means_long %>%
    left_join(obs_means, by = "decile") %>%
    group_by(decile) %>%
    summarise(pval = mean(value >= obs_value), .groups = "drop") %>%
    pull(pval)
  
  pval_sds <- ppc_sds_long %>%
    left_join(obs_sds, by = "decile") %>%
    group_by(decile) %>%
    summarise(pval = mean(value >= obs_value), .groups = "drop") %>%
    pull(pval)
  
  result <- matrix(
    c(pval_means, pval_sds),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      c("mean", "sd"),
      decile_levels
    )
  )
  
  result
}

#LC
lc_mean_list <- readRDS(paste0(input, "lc/lc_posterior_predictive_check_means.RDS"))
lc_sd_list <- readRDS(paste0(input, "lc/lc_posterior_predictive_check_sd.RDS"))

species<- c(
  "Homo_sapiens",
  "Neanderthal",
  "Pan_paniscus",
  "Pan_troglodytes",
  "Gorilla_beringei",
  "Gorilla_gorilla",
  "Pongo_abelii",
  "Pongo_pygmaeus"
)

total <- matrix(data = NA, nrow = 0, ncol = 9)
for (i in species) {
  p <- calc_ppc_pvalues(
    ppc_means    = lc_mean_list[[i]],
    ppc_sds      = lc_sd_list[[i]],
    observed_data = lc_data,
    species      = i
  )
  total <- rbind(total, cbind(c(i, i), p))
  # total <- rbind(total, p)
}
total <- as.data.frame(total)


for (i in species) {
  p <- plot_ppc(
    ppc_means    = lc_mean_list[[i]],
    ppc_sds      = lc_sd_list[[i]],
    observed_data = lc_data,
    species      = i
  )
  print(p)
  ggsave(
    paste0(output, "/ppc/", i, ".pdf"),
    plot   = p,
    width  = 14, height = 6,
    device = cairo_pdf
  )
}

### UI2
ui2_mean_list <- readRDS(paste0(input, "ui2/ui2_posterior_predictive_check_means.RDS"))
ui2_sd_list <- readRDS(paste0(input, "ui2/ui2_posterior_predictive_check_sd.RDS"))

species<- c(
  "Homo_sapiens",
  "Neanderthal",
  "Pan_paniscus",
  "Pan_troglodytes"
)

for (i in species) {
  p <- plot_ppc(
    ppc_means    = ui2_mean_list[[i]],
    ppc_sds      = ui2_sd_list[[i]],
    observed_data = ui2_data,
    species      = i
  )
  print(p)
  ggsave(
    paste0(output, "/ppc/", i, "_ui2.pdf"),
    plot   = p,
    width  = 14, height = 6,
    device = cairo_pdf
  )
}

total <- matrix(data = NA, nrow = 0, ncol = 9)
for (i in species) {
  p <- calc_ppc_pvalues(
    ppc_means    = ui2_mean_list[[i]],
    ppc_sds      = ui2_sd_list[[i]],
    observed_data = ui2_data,
    species      = i
  )
  total <- rbind(total, cbind(c(i, i), p))
}
total


# simulation study figure -------------------------------------------------

# Table #

parseTable <- function(cov_res_full, timing_res_full, cov_res_wo) {
  join_keys <- c("n_tips", "n_traits", "n_imp", "n_obs")
  
  cov_res_full %>%
    full_join(timing_res_full, by = join_keys) %>%
    full_join(cov_res_wo, by = join_keys, suffix = c("_full", "_wo")) %>%
    transmute(
      n_tips, n_traits, n_imp, n_obs,
      X    = mean_full,
      SE_X = se_full,
      Y    = mean_wo,
      SE_Y = se_wo,
      TIME = timing / 100
    )
}

full_model_coverage_result <- readRDS(paste0(input, "simulation_study/full_model_coverage_result_parsed.rds"))
full_model_timing_result <- readRDS(paste0(input, "simulation_study/full_model_timing_result_parsed.rds"))
wo_intra_coverage_result <- readRDS(paste0(input, "simulation_study/wo_intra_coverage_result_parsed.rds"))

table_data <- parseTable(full_model_coverage_result, full_model_timing_result, wo_intra_coverage_result)

CELL_FMT <- "\\cellentry{%s}{%s}{%s}{%s}{%s}"

fmt_val <- function(x, digits) {
  if (length(x) == 0 || is.na(x)) "--" else sprintf(paste0("%.", digits, "f"), x)
}

format_cell <- function(d, tips, traits, imp, obs) {
  row <- filter(d, n_tips == tips, n_traits == traits, n_imp == imp, n_obs == obs)
  
  if (nrow(row) == 0) {
    return(sprintf(CELL_FMT, "--", "--", "--", "--", "--"))
  }
  if (nrow(row) > 1) {
    stop(sprintf("Duplicate rows for n_tips=%s, n_traits=%s, n_imp=%s, n_obs=%s",
                 tips, traits, imp, obs))
  }
  
  full_ok <- !anyNA(row[c("X", "SE_X")])
  wo_ok   <- !anyNA(row[c("Y", "SE_Y")])
  
  sprintf(
    CELL_FMT,
    if (full_ok) fmt_val(row$X, 2)    else "--",
    if (full_ok) fmt_val(row$SE_X, 3) else "--",
    if (wo_ok)   fmt_val(row$Y, 2)    else "--",
    if (wo_ok)   fmt_val(row$SE_Y, 3) else "--",
    fmt_val(row$TIME, 1)
  )
}

obs_levels   <- c(2, 4, 8, 16)
imp_levels   <- c(0, 8, 16, 32)
trait_levels <- c(8, 16, 32)
tips_levels  <- c(8, 16, 32)

build_row <- function(d, tips, traits, imp) {
  cells <- vapply(obs_levels, function(o) format_cell(d, tips, traits, imp, o), character(1))
  sprintf("%d & %d & %s \\\\", traits, imp, paste(cells, collapse = " & "))
}

build_trait_block <- function(d, tips, traits) {
  rows <- vapply(imp_levels, function(i) build_row(d, tips, traits, i), character(1))
  paste(c(sprintf("%% %d traits", traits), rows), collapse = "\n")
}

build_panel <- function(d, tips, first = FALSE) {
  blocks <- vapply(trait_levels, function(tr) build_trait_block(d, tips, tr), character(1))
  body   <- paste(blocks, collapse = "\n\\midrule\n")
  
  comment <- paste(
    "% ------------------------------------------------------------------",
    sprintf("%% Panel: %d taxa", tips),
    "% ------------------------------------------------------------------",
    sep = "\n"
  )
  if (!first) comment <- paste(comment, "\\midrule", sep = "\n")  # separates from the previous panel
  
  header <- sprintf("\\multicolumn{6}{c}{\\textbf{%d Taxa}} \\\\\n\\midrule", tips)
  
  paste(comment, header, body, sep = "\n")
}

build_table_body <- function(d) {
  panels <- vapply(seq_along(tips_levels), function(i) {
    build_panel(d, tips_levels[i], first = (i == 1))
  }, character(1))
  paste(panels, collapse = "\n \n")
}

table_body <- build_table_body(table_data)

header_tex <- r"---(\clearpage
\begingroup
\scriptsize
\renewcommand{\arraystretch}{1.2}
\setlength{\tabcolsep}{3pt}
\begin{longtable}{cccccc}
\caption{Coverage probability and Markov chain Monte Carlo (MCMC) runtime from a full-factorial simulation study that varied the number of taxa (8, 16, or 32; shown top to bottom in separate panels), the number of traits (8, 16, or 32; panel rows), the number of observations per taxon (2, 4, 8, or 16; panel columns), and the number of missing observations across the entire simulated dataset (0, 8, 16, or 32; rows within each cell). For each combination of conditions we simulated 100 datasets and inferred parameters with both the full model, which accounts for intraspecific variation, and the taxon means model, which takes the mean value across observations in a taxon as known without uncertainty. Within each cell, the top line reports coverage probability (mean $\pm$ standard error across replicates) for the full model (left of the slash) and the taxon means model (right of the slash); the bottom line reports the mean wall clock time in minutes per ten million MCMC cycles for the full model running serially. The expected coverage probability is 95\%. } \label{table:coverage} \\
 
\toprule
 & & \multicolumn{4}{c}{Observations per Taxon} \\
\cmidrule(lr){3-6}
\textbf{Traits} & \textbf{Missing} & \textbf{2} & \textbf{4} & \textbf{8} & \textbf{16} \\
\midrule
\endfirsthead
 
\multicolumn{6}{l}{\textit{Table \thetable{} continued from previous page}} \\
\toprule
 & & \multicolumn{4}{c}{Observations per Taxon} \\
\cmidrule(lr){3-6}
\textbf{Traits} & \textbf{Missing} & \textbf{2} & \textbf{4} & \textbf{8} & \textbf{16} \\
\midrule
\endhead
 
\midrule
\multicolumn{6}{r}{\textit{continued on next page}} \\
\endfoot
 
\bottomrule
\endlastfoot
)---"

footer_tex <- r"---(
\end{longtable}
\endgroup
 
\bigskip
{\footnotesize\textit{Note.}~Each cell reports coverage probability (mean $\pm$ SE) for the full model / taxon-means model above the rule, and mean MCMC runtime (minutes per 10 million cycles, full model) below the rule.}
)---"

full_table_tex <- paste0(header_tex, "\n", table_body, "\n", footer_tex)

writeLines(full_table_tex, paste0(output, "/coverage_table.tex"))


# Runk and timing #
full_model_coverage_result <- readRDS(paste0(input, "simulation_study/full_model_coverage_result_parsed.rds"))
wo_intra_coverage_result   <- readRDS(paste0(input, "simulation_study/wo_intra_coverage_result_parsed.rds"))
full_model_timing_result   <- readRDS(paste0(input, "simulation_study/full_model_timing_result_parsed.rds"))

sim_join_keys <- c("n_tips", "n_traits", "n_imp", "n_obs")
sim_obs_levels <- c(2, 4, 8, 16)
sim_trait_levels <- c(8, 16, 32)
sim_tips_levels <- c(8, 16, 32)

sim_model_colors <- c(
  "Full model"        = "#1f78b4",
  "Taxon means model" = "#a31e22"
)

sim_trait_colors <- setNames(
  colorRampPalette(c("#f0b6b7", "#c14b4f", "#6d1417"))(length(sim_trait_levels)),
  paste(sim_trait_levels, "traits")
)

coverage_dat <- bind_rows(
  full_model_coverage_result %>%
    dplyr::select(all_of(sim_join_keys), coverage = mean) %>%
    mutate(model = "Full model"),
  wo_intra_coverage_result %>%
    dplyr::select(all_of(sim_join_keys), coverage = mean) %>%
    mutate(model = "Taxon means model")
) %>%
  drop_na(coverage) %>%
  mutate(
    coverage = coverage / if_else(max(coverage, na.rm = TRUE) > 1.5, 100, 1),
    model    = factor(model, levels = names(sim_model_colors)),
    obs_lab  = factor(n_obs, levels = sim_obs_levels)
  )

coverage_labels <- coverage_dat %>%
  group_by(obs_lab, model) %>%
  summarise(coverage = mean(coverage), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = coverage) %>%
  mutate(label = paste0(
    sprintf("%.1f", `Full model` * 100), "% / ",
    sprintf("%.1f", `Taxon means model` * 100), "%"
  ))

coveragePanel <- ggplot(coverage_dat, aes(x = obs_lab, y = coverage, fill = model)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  geom_half_violin(
    data = filter(coverage_dat, model == "Full model"),
    alpha = 0.7, scale = "width", side = "l"
  ) +
  geom_half_violin(
    data = filter(coverage_dat, model == "Taxon means model"),
    alpha = 0.7, scale = "width", side = "r"
  ) +
  geom_half_boxplot(
    data = filter(coverage_dat, model == "Full model"),
    alpha = 0.7, side = "l", width = 0.35, outlier.shape = NA
  ) +
  geom_half_boxplot(
    data = filter(coverage_dat, model == "Taxon means model"),
    alpha = 0.7, side = "r", width = 0.35, outlier.shape = NA
  ) +
  geom_text(
    data = coverage_labels,
    aes(x = obs_lab, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.5, size = 4, color = "grey30"
  ) +
  scale_fill_manual(values = sim_model_colors) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.14))
  ) +
  labs(
    x    = "Observations per taxon",
    y    = "Coverage probability",
    fill = "Model"
  ) +
  theme_minimal(base_family = "Georgia") +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 12)
  )

timing_dat <- full_model_timing_result %>%
  mutate(timing = timing / 100) %>%
  drop_na(timing) %>%
  group_by(n_tips, n_traits, n_obs) %>%
  summarise(
    mean_time = mean(timing),
    lo_time   = min(timing),
    hi_time   = max(timing),
    .groups   = "drop"
  ) %>%
  mutate(
    obs_lab   = factor(n_obs, levels = sim_obs_levels),
    trait_lab = factor(paste(n_traits, "traits"), levels = paste(sim_trait_levels, "traits")),
    tips_lab  = factor(paste(n_tips, "taxa"), levels = paste(sim_tips_levels, "taxa"))
  )

timingPanel <- ggplot(timing_dat, aes(x = obs_lab, y = mean_time, color = trait_lab, group = trait_lab)) +
  geom_linerange(aes(ymin = lo_time, ymax = hi_time), linewidth = 0.5, alpha = 0.8) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_wrap(~ tips_lab, nrow = 1) +
  scale_color_manual(values = sim_trait_colors) +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000),
    minor_breaks = NULL
  ) +
  labs(
    x     = "Observations per taxon",
    y     = "Mean wall clock time per 10 million\nMCMC cycles (minutes, log scale)",
    color = "Traits"
  ) +
  theme_minimal(base_family = "Georgia") +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(size = 12)
  )

sim_tag_theme <- theme(plot.tag = element_text(family = "Georgia", face = "bold", size = 13))

simulationFigure <- (coveragePanel + sim_tag_theme) / (timingPanel + sim_tag_theme) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")

simulationFigure
ggsave(paste0(output, "simulationStudy.svg"), plot = simulationFigure, width = 12, height = 9)

#### Table Latex Generation ####

tableTaxonOrder <- c(
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

tableDecileLabels <- paste0("Decile ", 3:10)

fmtTableNum <- function(x, digits = 2, big.mark = "") {
  vapply(x, function(v) {
    if (is.na(v)) {
      "-"
    } else {
      format(round(v, digits), big.mark = big.mark, trim = TRUE,
             scientific = FALSE, digits = 15)
    }
  }, character(1), USE.NAMES = FALSE)
}

orderByTaxon <- function(x, taxa = tableTaxonOrder) {
  setNames(as.numeric(x)[match(taxa, names(x))], taxa)
}

writeTableOutputs <- function(tex, dat, base) {
  writeLines(tex, paste0(output, base, ".tex"))
  saveRDS(dat, paste0(output, base, ".RDS"))
}

evoVcvMatrix <- function(fullPosterior, meansPosterior, cellFun) {
  fullCols  <- grep("^evo_vcv_", colnames(fullPosterior), value = TRUE)
  meansCols <- grep("^evo_vcv_", colnames(meansPosterior), value = TRUE)
  sharedCols <- fullCols[fullCols %in% meansCols]
  
  idx <- regmatches(sharedCols, regexpr("\\(\\d+,\\d+\\)", sharedCols))
  idx <- do.call(rbind, strsplit(gsub("[()]", "", idx), ","))
  rows <- as.integer(idx[, 1]) + 1
  cols <- as.integer(idx[, 2]) + 1
  
  out <- matrix(NA_real_, nrow = max(rows), ncol = max(cols),
                dimnames = list(tableDecileLabels, tableDecileLabels))
  for (k in seq_along(sharedCols)) {
    out[rows[k], cols[k]] <- cellFun(
      fullPosterior[[sharedCols[k]]],
      meansPosterior[[sharedCols[k]]]
    )
  }
  out
}

matrixBlockTex <- function(mat, header) {
  rows <- vapply(seq_len(nrow(mat)), function(i) {
    paste0("    ", rownames(mat)[i], " & ",
           paste(fmtTableNum(mat[i, ]), collapse = " & "), " \\\\")
  }, character(1))
  
  paste(
    c(
      paste0("    \\multicolumn{1}{l}{", header, "} & \\\\"),
      "    \\midrule",
      paste0("    & ", paste(colnames(mat), collapse = " & "), " \\\\"),
      "    \\midrule",
      rows
    ),
    collapse = "\n"
  )
}

calcKLDivergenceInverseWishart <- function(scalePost, dofPost, scalePrior, dofPrior) {
  p <- 8
  V1 <- solve(scalePost)
  V2 <- solve(scalePrior)
  n1 <- dofPost
  n2 <- dofPrior
  term1 <- n2 * as.numeric(
    determinant(V2, logarithm = TRUE)$modulus -
      determinant(V1, logarithm = TRUE)$modulus
  )
  term2 <- n1 * sum(diag(solve(V2) %*% V1))
  term3 <- 2 * (CholWishart::lmvgamma(n2 / 2, p) - CholWishart::lmvgamma(n1 / 2, p))
  term4 <- (n1 - n2) * CholWishart::mvdigamma(n1 / 2, p)
  term5 <- -n1 * p
  
  as.numeric(0.5 * (term1 + term2 + term3 + term4 + term5))
}

calcSymmetrizedKLDivergence <- function(posteriorFit1, posteriorFit2) {
  klforward <- calcKLDivergenceInverseWishart(
    scalePost  = posteriorFit1$scale,
    dofPost    = posteriorFit1$nu,
    scalePrior = posteriorFit2$scale,
    dofPrior   = posteriorFit2$nu
  )
  klbackward <- calcKLDivergenceInverseWishart(
    scalePost  = posteriorFit2$scale,
    dofPost    = posteriorFit2$nu,
    scalePrior = posteriorFit1$scale,
    dofPrior   = posteriorFit1$nu
  )
  
  round(klforward + klbackward, 2)
}

## table:klDivergences

lc_posteriorFits  <- readRDS(paste0(input, "lc/lc_dec3_10_posterior_fits.RDS"))
ui2_posteriorFits <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_posterior_fits.RDS"))

klTraits    <- 8
klPriorDOF  <- 10
klPriorScale <- matrix(1e-6, klTraits, klTraits)
diag(klPriorScale) <- 1.0

klPriorToPosterior <- function(postFits, analysis) {
  kl <- vapply(postFits, function(fit) {
    calcKLDivergenceInverseWishart(
      scalePost  = fit$scale,
      dofPost    = fit$nu,
      scalePrior = klPriorScale,
      dofPrior   = klPriorDOF
    )
  }, numeric(1))
  nu <- vapply(postFits, function(fit) as.numeric(fit$nu), numeric(1))
  
  data.frame(
    analysis = analysis,
    taxon    = tableTaxonOrder,
    kl       = unname(orderByTaxon(kl)),
    nu       = unname(orderByTaxon(nu)),
    stringsAsFactors = FALSE
  )
}

klDivergencesData <- rbind(
  klPriorToPosterior(lc_posteriorFits, "C1"),
  klPriorToPosterior(ui2_posteriorFits, "I2")
)

klDivergencesRow <- function(dat, analysis, label) {
  vals  <- dat[dat$analysis == analysis, ]
  vals  <- vals[match(tableTaxonOrder, vals$taxon), ]
  cells <- paste0(fmtTableNum(vals$kl), " (", fmtTableNum(vals$nu), ")")
  cells[is.na(vals$kl)] <- "-"
  paste0("    ", label, " & ", paste(cells, collapse = " & "), " \\\\")
}

klDivergencesTex <- paste(
  r"---(\begin{table*}[p]
    \centering
    \resizebox{\textwidth}{!}{
    \setlength{\tabcolsep}{2pt}
    \begin{tabular}{lccccccccc}
    \toprule
          \rotatebox{45}{Analysis} 
        & \rotatebox{45}{Evolutionary VCV}
        & \rotatebox{45}{Modern human VCV}
        & \rotatebox{45}{Neandertal VCV}
        & \rotatebox{45}{\textit{Pan paniscus} VCV}
        & \rotatebox{45}{\textit{Pan troglodytes} VCV}
        & \rotatebox{45}{\textit{G. beringei} VCV}
        & \rotatebox{45}{\textit{G. gorilla} VCV}
        & \rotatebox{45}{\textit{Po. abelii} VCV}
        & \rotatebox{45}{\textit{Po. pygmaeus} VCV} \\
    \midrule)---",
  klDivergencesRow(klDivergencesData, "C1", "C\\textsubscript{1}"),
  "    \\midrule",
  klDivergencesRow(klDivergencesData, "I2", "I\\textsuperscript{2}"),
  r"---(    \bottomrule
    \end{tabular}}
    \caption{Kullback-Leibler divergences from the prior to the posterior distribution for the variance-covariance matrices. These posterior distributions are Inverse Wishart distributed; the maximum likelihood estimate for each posterior distribution's degree of freedom is listed in parentheses. Lower degrees of freedom indicate heavier tails.}
    \label{table:klDivergences}
\end{table*})---",
sep = "\n"
)

writeTableOutputs(klDivergencesTex, klDivergencesData, "table_klDivergences")

## table:symklDivergences

lc_posteriorFits  <- readRDS(paste0(input, "lc/lc_dec3_10_posterior_fits.RDS"))
ui2_posteriorFits <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_posterior_fits.RDS"))

symklComparisons <- list(
  list(analysis = "C1", label = "Modern human vs. Neandertals",
       taxa = c("Homo_sapiens", "Neanderthal")),
  list(analysis = "C1", label = "\\textit{Pan paniscus} vs. \\textit{Pan troglodytes}",
       taxa = c("Pan_troglodytes", "Pan_paniscus")),
  list(analysis = "C1", label = "\\textit{Gorilla beringei} vs. \\textit{Gorilla gorilla}",
       taxa = c("Gorilla_beringei", "Gorilla_gorilla")),
  list(analysis = "C1", label = "\\textit{Pongo abelii} vs. \\textit{Pongo pygmaeus}",
       taxa = c("Pongo_abelii", "Pongo_pygmaeus")),
  list(analysis = "I2", label = "Modern human vs. Neandertals",
       taxa = c("Homo_sapiens", "Neanderthal")),
  list(analysis = "I2", label = "\\textit{Pan paniscus} vs. \\textit{Pan troglodytes}",
       taxa = c("Pan_troglodytes", "Pan_paniscus"))
)

symklFits <- list(C1 = lc_posteriorFits, I2 = ui2_posteriorFits)

symklDivergencesData <- do.call(rbind, lapply(symklComparisons, function(cmp) {
  fits <- symklFits[[cmp$analysis]]
  data.frame(
    analysis   = cmp$analysis,
    comparison = cmp$label,
    taxon1     = cmp$taxa[1],
    taxon2     = cmp$taxa[2],
    symkl      = calcSymmetrizedKLDivergence(fits[[cmp$taxa[1]]], fits[[cmp$taxa[2]]]),
    stringsAsFactors = FALSE
  )
}))

symklDivergencesBlock <- function(dat, analysis, header) {
  vals <- dat[dat$analysis == analysis, ]
  rows <- vapply(seq_len(nrow(vals)), function(i) {
    paste0("    ", vals$comparison[i], "  & ", fmtTableNum(vals$symkl[i]), " \\\\\n    \\midrule")
  }, character(1))
  paste(c(paste0("    \\multicolumn{1}{l}{", header, "} & \\\\"), "    \\midrule", rows),
        collapse = "\n")
}

symklDivergencesTex <- paste(
  r"---(\begin{table*}[p]
    \centering
    \resizebox{\textwidth}{!}{
    \begin{tabular}{lc}
    \toprule
    \multicolumn{1}{l}{Comparison} & Symmetrized KL divergence \\
    \midrule)---",
  symklDivergencesBlock(symklDivergencesData, "C1", "\\textit{C\\textsubscript{1}}"),
  symklDivergencesBlock(symklDivergencesData, "I2", "\\textit{I\\textsuperscript{2}}"),
  r"---(    \bottomrule
    \end{tabular}}
    \caption{Symmetrized Kullback-Leibler divergences measuring the divergence in the inferred C\textsubscript{1} and I\textsuperscript{2} intraspecific variance-covariance distributions between members of a genus.}
    \label{table:symklDivergences}
\end{table*})---",
sep = "\n"
)

writeTableOutputs(symklDivergencesTex, symklDivergencesData, "table_symklDivergences")

## table:frechetVariances
lc_Frechet_Var  <- readRDS(paste0(input, "lc/lc_vcv_frechet_var.RDS"))
ui2_Frechet_Var <- readRDS(paste0(input, "ui2/ui2_vcv_frechet_var.RDS"))

lc_Frechet_Var_TaxonMeans  <- readRDS(paste0(input, "lc/lc_evo_vcv_frechet_var_species_mean.RDS"))
ui2_Frechet_Var_TaxonMeans <- readRDS(paste0(input, "ui2/ui2_evo_vcv_frechet_var_species_mean.RDS"))

# tableTaxonOrder is assumed to start with "evolutionary" followed by the 8 taxa,
# as used previously; strip "evolutionary" out to get just the intraspecific taxa
intraTaxonOrder <- tableTaxonOrder[tableTaxonOrder != "evolutionary"]

frechetVariancesData <- rbind(
  data.frame(analysis = "C1", taxon = tableTaxonOrder,
             frechet = unname(orderByTaxon(lc_Frechet_Var)), stringsAsFactors = FALSE),
  data.frame(analysis = "I2", taxon = tableTaxonOrder,
             frechet = unname(orderByTaxon(ui2_Frechet_Var)), stringsAsFactors = FALSE),
  data.frame(analysis = c("C1", "I2"),
             taxon = "evolutionary_taxonMeans",
             frechet = c(unname(lc_Frechet_Var_TaxonMeans), unname(ui2_Frechet_Var_TaxonMeans)),
             stringsAsFactors = FALSE)
)

frechetVariancesRow <- function(dat, analysis, label) {
  fullModel  <- dat[dat$analysis == analysis & dat$taxon == "evolutionary", "frechet"]
  taxonMeans <- dat[dat$analysis == analysis & dat$taxon == "evolutionary_taxonMeans", "frechet"]
  
  intra <- dat[dat$analysis == analysis & dat$taxon %in% intraTaxonOrder, ]
  intra <- intra[match(intraTaxonOrder, intra$taxon), ]
  
  vals <- c(fullModel, taxonMeans, intra$frechet)
  paste0("    ", label, " & ", paste(fmtTableNum(vals), collapse = " & "), " \\\\")
}

frechetVariancesTex <- paste(
  r"---(\begin{table*}[t]
    \centering
    \resizebox{\textwidth}{!}{
    \setlength{\tabcolsep}{4pt}
    \begin{tabular}{l cc cccccccc}
    \toprule
     & \multicolumn{2}{c}{Evolutionary VCV} & \multicolumn{8}{c}{Intraspecific VCV (full model)} \\
    \cmidrule(lr){2-3} \cmidrule(lr){4-11}
    \rotatebox{45}{Analysis}
        & \rotatebox{45}{Full model}
        & \rotatebox{45}{Taxon means model}
        & \rotatebox{45}{Modern human}
        & \rotatebox{45}{Neandertal}
        & \rotatebox{45}{\textit{Pan paniscus}}
        & \rotatebox{45}{\textit{Pan troglodytes}}
        & \rotatebox{45}{\textit{G. beringei}}
        & \rotatebox{45}{\textit{G. gorilla}}
        & \rotatebox{45}{\textit{Po. abelii}}
        & \rotatebox{45}{\textit{Po. pygmaeus}} \\
    \midrule)---",
  frechetVariancesRow(frechetVariancesData, "C1", "C\\textsubscript{1}"),
  "    \\midrule",
  frechetVariancesRow(frechetVariancesData, "I2", "I\\textsuperscript{2}"),
  r"---(    \bottomrule
    \end{tabular}}
    \caption{Fréchet variance for each VCV posterior distribution. The evolutionary VCV is compared
    between the full model, which accounts for intraspecific variation, and the taxon means model, in
    which species means are fixed at their observed sample means and only the evolutionary VCV is
    estimated; the remaining columns give the intraspecific VCV inferred for each taxon under the full
    model, which has no analog under the taxon means model. Higher values indicate greater uncertainty
    in the inferred VCV matrix.}
    \label{table:frechetVariances}
\end{table*})---",
sep = "\n"
)
writeTableOutputs(frechetVariancesTex, frechetVariancesData, "table_frechetVariances")

## table:PostPred

lc_data  <- read.csv("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv")
ui2_data <- read.csv("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/UI2dec3_10_no_pongo.csv")

postPredFiles <- c(
  Homo_sapiens     = "hsPostPred.rds",
  Neanderthal      = "neanderthalPostPred.rds",
  Pan_paniscus     = "panpaniscusPostPred.rds",
  Pan_troglodytes  = "pantroglodytesPostPred.rds",
  Gorilla_beringei = "gorrillaberingeiPostPred.rds",
  Gorilla_gorilla  = "gorillagorillaPostPred.rds",
  Pongo_abelii     = "pongoabeliiPostPred.rds",
  Pongo_pygmaeus   = "pongopygmaeusPostPred.rds"
)

postPredGroups <- list(
  list(label = "Modern humans",              taxon = "Homo_sapiens",     teeth = c("C1", "I2")),
  list(label = "Neandertals",                taxon = "Neanderthal",      teeth = c("C1", "I2")),
  list(label = "\\textit{Pan paniscus}",     taxon = "Pan_paniscus",     teeth = c("C1", "I2")),
  list(label = "\\textit{Pan troglodytes}",  taxon = "Pan_troglodytes",  teeth = c("C1", "I2")),
  list(label = "\\textit{Gorilla beringei}", taxon = "Gorilla_beringei", teeth = "C1"),
  list(label = "\\textit{Gorilla gorilla}",  taxon = "Gorilla_gorilla",  teeth = "C1"),
  list(label = "\\textit{Pongo abelii}",     taxon = "Pongo_abelii",     teeth = "C1"),
  list(label = "\\textit{Pongo pygmaeus}",   taxon = "Pongo_pygmaeus",   teeth = "C1")
)

postPredDirs   <- c(C1 = "lc/posteriorPredictive/", I2 = "ui2/posteriorPredictive/")
postPredCounts <- list(C1 = lc_data, I2 = ui2_data)

postPredSummary <- function(taxon, tooth) {
  pred <- readRDS(paste0(input, postPredDirs[[tooth]], postPredFiles[[taxon]]))
  data.frame(
    taxon  = taxon,
    tooth  = tooth,
    n      = sum(postPredCounts[[tooth]]$genus == taxon),
    decile = tableDecileLabels,
    mean   = vapply(tableDecileLabels, function(d) mean(pred[[d]]), numeric(1)),
    var    = vapply(tableDecileLabels, function(d) var(pred[[d]]), numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

postPredData <- do.call(rbind, lapply(postPredGroups, function(grp) {
  do.call(rbind, lapply(grp$teeth, function(tooth) postPredSummary(grp$taxon, tooth)))
}))

postPredToothLabels <- c(C1 = "C\\textsubscript{1}", I2 = "I\\textsuperscript{2}")

postPredCells <- function(dat, taxon, tooth) {
  vals <- dat[dat$taxon == taxon & dat$tooth == tooth, ]
  vals <- vals[match(tableDecileLabels, vals$decile), ]
  paste(paste0(fmtTableNum(vals$mean), " (", fmtTableNum(vals$var), ")"), collapse = " & ")
}

postPredRows <- function(dat, grp) {
  if (length(grp$teeth) == 1) {
    vals <- dat[dat$taxon == grp$taxon & dat$tooth == grp$teeth, ]
    return(paste0("    ", grp$label, " & ", vals$n[1], " & ",
                  postPredCells(dat, grp$taxon, grp$teeth), " \\\\"))
  }
  
  subRows <- vapply(grp$teeth, function(tooth) {
    vals <- dat[dat$taxon == grp$taxon & dat$tooth == tooth, ]
    paste0("    \\multicolumn{1}{r}{", postPredToothLabels[[tooth]], "} & ", vals$n[1],
           " & ", postPredCells(dat, grp$taxon, tooth), " \\\\")
  }, character(1))
  
  paste(c(paste0("    ", grp$label, " & \\\\"), subRows), collapse = "\n")
}

postPredBody <- paste(
  c(
    postPredRows(postPredData, postPredGroups[[1]]),
    postPredRows(postPredData, postPredGroups[[2]]),
    "    \\midrule",
    postPredRows(postPredData, postPredGroups[[3]]),
    postPredRows(postPredData, postPredGroups[[4]]),
    "    \\midrule",
    postPredRows(postPredData, postPredGroups[[5]]),
    postPredRows(postPredData, postPredGroups[[6]]),
    "    \\midrule",
    postPredRows(postPredData, postPredGroups[[7]]),
    postPredRows(postPredData, postPredGroups[[8]])
  ),
  collapse = "\n"
)

postPredTex <- paste(
  paste0(
    r"---(\begin{table*}[h]
    \centering
    \resizebox{\textwidth}{!}{
    \begin{tabular}{lccccccccc}
    \toprule
    Group & Sample size & )---",
    paste(tableDecileLabels, collapse = " & "),
    " \\\\\n    \\midrule"
  ),
  postPredBody,
  r"---(    \bottomrule
    \end{tabular}}
    \caption{Comparison of posterior predictive perikymata per millimeter per decile means and variances (in parentheses). For modern humans, Neandertals, and \textit{Pan}, where we have both C\textsubscript{1} and I\textsuperscript{2} data, both values are listed. We only have C\textsubscript{1} data for \textit{Gorilla} and \textit{Pongo}.}
    \label{table:PostPred}
\end{table*})---",
sep = "\n"
)

writeTableOutputs(postPredTex, postPredData, "table_PostPred")

## table:overlapEmpiricalFull

lc_posterior <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10.tsv")))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ]

lc_posterior_species_means <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_species_means.tsv")))
lc_posterior_species_means <- lc_posterior_species_means[round(0.1 * nrow(lc_posterior_species_means)) : nrow(lc_posterior_species_means), ]

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ]

ui2_posterior_species_means <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_species_means.tsv")))
ui2_posterior_species_means <- ui2_posterior_species_means[round(0.1 * nrow(ui2_posterior_species_means)) : nrow(ui2_posterior_species_means), ]

evoVcvOverlap <- function(x, y) {
  100 * overlapping::overlap(list(x, y))$OV[[1]]
}

overlapEmpiricalFullData <- list(
  C1 = evoVcvMatrix(lc_posterior, lc_posterior_species_means, evoVcvOverlap),
  I2 = evoVcvMatrix(ui2_posterior, ui2_posterior_species_means, evoVcvOverlap)
)

overlapEmpiricalFullTex <- paste(
  r"---(\begin{table*}[p]
    \centering
    \resizebox{\textwidth}{!}{
    \begin{tabular}{lcccccccc}
    \toprule)---",
  matrixBlockTex(overlapEmpiricalFullData$C1, "\\textit{C\\textsubscript{1}}"),
  "    \\midrule",
  matrixBlockTex(overlapEmpiricalFullData$I2, "\\textit{I\\textsuperscript{2}}"),
  r"---(    \bottomrule
    \end{tabular}
    }
    \caption{Overlap percentage for each element of the evolutionary variance-covariance matrices between the posterior distributions inferred with the full model and the empirical species means model. Higher overlap indicates the posterior probability distribution for that element of the evolutionary VCV matrix is more similar between the two models.}
    \label{table:overlapEmpiricalFull}
\end{table*})---",
sep = "\n"
)

writeTableOutputs(overlapEmpiricalFullTex, overlapEmpiricalFullData, "table_overlapEmpiricalFull")

## table:fullEmpiricalVar

lc_posterior <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10.tsv")))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ]

lc_posterior_species_means <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_species_means.tsv")))
lc_posterior_species_means <- lc_posterior_species_means[round(0.1 * nrow(lc_posterior_species_means)) : nrow(lc_posterior_species_means), ]

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ]

ui2_posterior_species_means <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_species_means.tsv")))
ui2_posterior_species_means <- ui2_posterior_species_means[round(0.1 * nrow(ui2_posterior_species_means)) : nrow(ui2_posterior_species_means), ]

evoVcvMapDifference <- function(x, y) {
  bayestestR::map_estimate(x)$MAP_Estimate - bayestestR::map_estimate(y)$MAP_Estimate
}

fullEmpiricalVarData <- list(
  C1 = evoVcvMatrix(lc_posterior, lc_posterior_species_means, evoVcvMapDifference),
  I2 = evoVcvMatrix(ui2_posterior, ui2_posterior_species_means, evoVcvMapDifference)
)

fullEmpiricalVarTex <- paste(
  r"---(\begin{table*}[p]
    \centering
    \resizebox{\textwidth}{!}{
    \begin{tabular}{lcccccccc}
    \toprule)---",
  matrixBlockTex(fullEmpiricalVarData$C1, "\\textit{C\\textsubscript{1}}"),
  "    \\midrule",
  matrixBlockTex(fullEmpiricalVarData$I2, "\\textit{I\\textsuperscript{2}}"),
  r"---(    \bottomrule
    \end{tabular}
    }
    \caption{Difference in maximum \textit{a posteriori} estimate for each element of the evolutionary variance-covariance matrices in the posterior distributions inferred with the full model and the empirical species means model. Positive values indicate the full model had a greater maximum \textit{a posteriori} estimate for that element.}
    \label{table:fullEmpiricalVar}
\end{table*})---",
sep = "\n"
)

writeTableOutputs(fullEmpiricalVarTex, fullEmpiricalVarData, "table_fullEmpiricalVar")

## table:meanVariances

lc_posterior <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10.tsv")))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ]

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ]

calcPosteriorMeanVariance <- function(posterior_df) {
  taxa <- tableTaxonOrder[-1]
  
  present <- taxa[vapply(taxa, function(taxon) {
    any(str_detect(names(posterior_df), paste0("^", taxon, "_mean_\\d+$")))
  }, logical(1))]
  
  vapply(present, function(taxon) {
    meanCols <- names(posterior_df)[str_detect(names(posterior_df),
                                               paste0("^", taxon, "_mean_\\d+$"))]
    meanCols <- meanCols[order(as.integer(str_extract(meanCols, "\\d+$")))]
    sum(sapply(posterior_df[, meanCols], var))
  }, numeric(1))
}

meanVariancesData <- rbind(
  data.frame(analysis = "C1", taxon = tableTaxonOrder[-1],
             meanVar = unname(orderByTaxon(calcPosteriorMeanVariance(lc_posterior),
                                           tableTaxonOrder[-1])),
             stringsAsFactors = FALSE),

  data.frame(analysis = "I2", taxon = tableTaxonOrder[-1],
             meanVar = unname(orderByTaxon(calcPosteriorMeanVariance(ui2_posterior),
                                           tableTaxonOrder[-1])),
             stringsAsFactors = FALSE)
)

meanVariancesRow <- function(dat, analysis, label) {
  vals <- dat[dat$analysis == analysis, ]
  vals <- vals[match(tableTaxonOrder[-1], vals$taxon), ]
  paste0("    ", label, " & ", paste(fmtTableNum(vals$meanVar), collapse = " & "), " \\\\")
}

meanVariancesTex <- paste(
  r"---(\begin{table*}[p]
    \centering
    \resizebox{\textwidth}{!}{
    \setlength{\tabcolsep}{1pt}
    \begin{tabular}{lcccccccc}
    \toprule
              \rotatebox{45}{Analysis}
            & \rotatebox{45}{Modern human $\sigma^2_\mu$}
            & \rotatebox{45}{Neandertal $\sigma^2_\mu$}
            & \rotatebox{45}{\textit{Pan paniscus} $\sigma^2_\mu$}
            & \rotatebox{45}{\textit{Pan troglodytes} $\sigma^2_\mu$}
            & \rotatebox{45}{\textit{G. beringei} $\sigma^2_\mu$}
            & \rotatebox{45}{\textit{G. gorilla} $\sigma^2_\mu$}
            & \rotatebox{45}{\textit{Po. abelii} $\sigma^2_\mu$}
            & \rotatebox{45}{\textit{Po. pygmaeus} $\sigma^2_\mu$} \\
    \midrule)---",
  meanVariancesRow(meanVariancesData, "C1", "C\\textsubscript{1}"),
  "    \\midrule",
  meanVariancesRow(meanVariancesData, "I2", "I\\textsuperscript{2}"),
  r"---(    \bottomrule
    \end{tabular}}
    \caption{Total variance in the posterior distribution for the inferred intraspecific means ($\sigma^2_\mu$). Calculated as the sum across all deciles of the variance in the posterior distribution for that decile.}
    \label{table:meanVariances}
\end{table*})---",
sep = "\n"
)

writeTableOutputs(meanVariancesTex, meanVariancesData, "table_meanVariances")

## table:symKLevoIntra

lc_posteriorFits             <- readRDS(paste0(input, "lc/lc_dec3_10_posterior_fits.RDS"))
lc_no_hominin_posteriorFits  <- readRDS(paste0(input, "lc/lc_dec3_10_no_hominin_posterior_fits.RDS"))
ui2_posteriorFits            <- readRDS(paste0(input, "ui2/ui2_dec3_10_no_pongo_posterior_fits.RDS"))

symKLevoIntraValues <- function(postFits) {
  taxa <- tableTaxonOrder[-1]
  taxa <- taxa[taxa %in% names(postFits)]
  vapply(taxa, function(taxon) {
    calcSymmetrizedKLDivergence(postFits$evolutionary, postFits[[taxon]])
  }, numeric(1))
}

symKLevoIntraData <- rbind(
  data.frame(analysis = "C1", taxon = tableTaxonOrder[-1],
             symkl = unname(orderByTaxon(symKLevoIntraValues(lc_posteriorFits),
                                         tableTaxonOrder[-1])),
             stringsAsFactors = FALSE),
  data.frame(analysis = "C1_no_hominin", taxon = tableTaxonOrder[-1],
             symkl = unname(orderByTaxon(symKLevoIntraValues(lc_no_hominin_posteriorFits),
                                         tableTaxonOrder[-1])),
             stringsAsFactors = FALSE),
  data.frame(analysis = "I2", taxon = tableTaxonOrder[-1],
             symkl = unname(orderByTaxon(symKLevoIntraValues(ui2_posteriorFits),
                                         tableTaxonOrder[-1])),
             stringsAsFactors = FALSE)
)

symKLevoIntraRow <- function(dat, analysis, label) {
  vals <- dat[dat$analysis == analysis, ]
  vals <- vals[match(tableTaxonOrder[-1], vals$taxon), ]
  paste0("    ", label, " & ",
         paste(fmtTableNum(vals$symkl, big.mark = ","), collapse = " & "), " \\\\")
}

symKLevoIntraTex <- paste(
  r"---(\begin{table*}[p]
    \centering
    \resizebox{\textwidth}{!}{
    \setlength{\tabcolsep}{1pt}
    \begin{tabular}{lcccccccc}
    \toprule
              \rotatebox{45}{Analysis}
            & \rotatebox{45}{Modern human VCV}
            & \rotatebox{45}{Neandertal VCV}
            & \rotatebox{45}{\textit{Pan paniscus} VCV}
            & \rotatebox{45}{\textit{Pan troglodytes} VCV}
            & \rotatebox{45}{\textit{G. beringei} VCV}
            & \rotatebox{45}{\textit{G. gorilla} VCV}
            & \rotatebox{45}{\textit{Po. abelii} VCV}
            & \rotatebox{45}{\textit{Po. pygmaeus} VCV} \\
    \midrule)---",
  symKLevoIntraRow(symKLevoIntraData, "C1", "C\\textsubscript{1}"),
  "    \\midrule",
  symKLevoIntraRow(symKLevoIntraData, "C1_no_hominin", "C\\textsubscript{1} w/o h"),
  "    \\midrule",
  symKLevoIntraRow(symKLevoIntraData, "I2", "I\\textsuperscript{2}"),
  r"---(    \bottomrule
    \end{tabular}}
    \caption{Symmetrized Kullback-Leibler divergences measuring the divergence in the inferred evolutionary VCV matrix with each intraspecific VCV matrix. }
    \label{table:symKLevoIntra}
\end{table*})---",
sep = "\n"
)

writeTableOutputs(symKLevoIntraTex, symKLevoIntraData, "table_symKLevoIntra")
