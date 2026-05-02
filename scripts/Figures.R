library(ape)
library(bayestestR)
library(data.table)
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

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::mutate)

input <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/"
output <- "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/figures/"

lc_posterior <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10.tsv")))
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

lc_posterior_no_hominin <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_no_hominin.tsv")))
lc_posterior_no_hominin  <- lc_posterior_no_hominin[round(0.1 * nrow(lc_posterior_no_hominin)) : nrow(lc_posterior_no_hominin), ] #apply burnin

lc_posterior_species_means <- as.data.frame(fread(paste0(input, "lc/lc_dec3_10_species_means.tsv")))
lc_posterior_species_means  <- lc_posterior_species_means[round(0.1 * nrow(lc_posterior_species_means)) : nrow(lc_posterior_species_means), ] #apply burnin

ui2_posterior <- as.data.frame(fread(paste0(input, "ui2/ui2_dec3_10_no_pongo.tsv")))
ui2_posterior <- ui2_posterior[round(0.1 * nrow(ui2_posterior)) : nrow(ui2_posterior), ] #apply burnin

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
map_estimate <- function(x) {
  as.numeric(bayestestR::map_estimate(x))
}

#### LC with hominins ####
#evo VCV MAP
evo_vcv_cols <- paste0("evo_vcv_(", rep(0:7, each = 8), ",", rep(0:7, times = 8), ")")

shared_max <- max(
  dplyr::select(lc_posterior, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist(),
  dplyr::select(lc_posterior_no_hominin, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist()
)

shared_min <- min(
  dplyr::select(lc_posterior, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist(),
  dplyr::select(lc_posterior_no_hominin, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist()
  )

evo_map <- dplyr::select(lc_posterior, all_of(evo_vcv_cols)) |>
  summarise(across(everything(), map_estimate)) |>
  pivot_longer(everything(), names_to = "element", values_to = "map") |>
  mutate(
    row = as.integer(sub(".*\\((\\d+),(\\d+)\\)", "\\1", element)),
    col = as.integer(sub(".*\\((\\d+),(\\d+)\\)", "\\2", element))
  )

decile_labels <- paste0("Decile ", 3:10)

evo_map <- evo_map |>
  mutate(
    row_label = factor(decile_labels[row + 1], levels = decile_labels),
    col_label = factor(decile_labels[col + 1], levels = decile_labels)
  )

evo_map$is_diag <- evo_map$col_label == evo_map$row_label

evoVCV <- ggplot(evo_map, aes(x = col_label, y = fct_rev(row_label), fill = map)) +
  geom_tile() +
  geom_tile(data = subset(evo_map, is_diag), color = "black", linewidth = 1.5, fill = NA) +
  geom_text(aes(label = round(map, 2)), size = 3, color = "black") +
  scale_fill_gradient(
    low    = "white",
    high   = "#a31e22",
    limits = c(shared_min, shared_max)
  )+
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
treeplot <- ggtree::rotate(treeplot, 14)
treeplot

#heatmap at tips
tip_order <- treeplot$data |>
  dplyr::filter(isTip) |>
  dplyr::arrange(y) |>
  dplyr::pull(label)

species_map <-c(
  "Pongo abelii" = "Pongo_abelii", 
  "Pongo pygmaeus" = "Pongo_pygmaeus", 
  "Gorilla beringei" = "Gorilla_beringei", 
  "Gorilla gorilla" =  "Gorilla_gorilla", 
  "Pan troglodytes" = "Pan_troglodytes", 
  "Pan paniscus" = "Pan_paniscus", 
  "Neanderthal" ="Neanderthal", 
  "Modern humans" = "Homo_sapiens"
)

mean_map <- lapply(names(species_map), function(tip_label) {
  sp        <- unname(species_map[tip_label])
  mean_cols <- paste0(sp, "_mean_", 0:7)
  vals      <- dplyr::select(lc_posterior, all_of(mean_cols)) |>
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


p1 <-  evoVCV +  treeplot+ tipplot
p1

#### LC without hominins ####
#evo VCV MAP
evo_map <- dplyr::select(lc_posterior_no_hominin, all_of(evo_vcv_cols)) |>
  summarise(across(everything(), map_estimate)) |>
  pivot_longer(everything(), names_to = "element", values_to = "map") |>
  mutate(
    row = as.integer(sub(".*\\((\\d+),(\\d+)\\)", "\\1", element)),
    col = as.integer(sub(".*\\((\\d+),(\\d+)\\)", "\\2", element))
  )

decile_labels <- paste0("Decile ", 3:10)

evo_map <- evo_map |>
  mutate(
    row_label = factor(decile_labels[row + 1], levels = decile_labels),
    col_label = factor(decile_labels[col + 1], levels = decile_labels)
  )

evo_map$is_diag <- evo_map$col_label == evo_map$row_label

evoVCV <- ggplot(evo_map, aes(x = col_label, y = fct_rev(row_label), fill = map)) +
  geom_tile() +
  geom_tile(data = subset(evo_map, is_diag), color = "black", linewidth = 1.5, fill = NA) +
  geom_text(aes(label = round(map, 2)), size = 3, color = "black") +
  scale_fill_gradient(
    low    = "white",
    high   = "#a31e22",
    limits = c(shared_min, shared_max)
  )+
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
plottree <- drop.tip(plottree, "Neanderthal")
plottree <- drop.tip(plottree, "Homo_sapiens")
plottree$tip.label <- gsub("_", " ", plottree$tip.label)
treeplot <- ggtree(plottree) + 
  geom_tiplab(aes(fontface = ifelse(label %in% c("Modern humans", "Neanderthal"), 2, 4)), family = "Georgia") +
  hexpand(0.55)
treeplot

#heatmap at tips
tip_order <- treeplot$data |>
  dplyr::filter(isTip) |>
  dplyr::arrange(y) |>
  dplyr::pull(label)

species_map <-c(
  "Pongo abelii" = "Pongo_abelii", 
  "Pongo pygmaeus" = "Pongo_pygmaeus", 
  "Gorilla beringei" = "Gorilla_beringei", 
  "Gorilla gorilla" =  "Gorilla_gorilla", 
  "Pan troglodytes" = "Pan_troglodytes", 
  "Pan paniscus" = "Pan_paniscus"
)

mean_map <- lapply(names(species_map), function(tip_label) {
  sp        <- unname(species_map[tip_label])
  mean_cols <- paste0(sp, "_mean_", 0:7)
  vals      <- dplyr::select(lc_posterior_no_hominin, all_of(mean_cols)) |>
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


p2 <-  evoVCV +treeplot + tipplot
p2

p3 <- p1/p2
p3

ggsave(paste0(output, "treeMAP.svg"), plot = p3, width = 14, height = 8)

# UI2 evo VCV ---------------------------------------------------------------

evo_vcv_cols <- paste0("evo_vcv_(", rep(0:7, each = 8), ",", rep(0:7, times = 8), ")")

shared_max <- max(
  dplyr::select(ui2_posterior, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist())

shared_min <- min(
  dplyr::select(ui2_posterior, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist()
)

evo_map <- dplyr::select(ui2_posterior, all_of(evo_vcv_cols)) |>
  summarise(across(everything(), map_estimate)) |>
  pivot_longer(everything(), names_to = "element", values_to = "map") |>
  mutate(
    row = as.integer(sub(".*\\((\\d+),(\\d+)\\)", "\\1", element)),
    col = as.integer(sub(".*\\((\\d+),(\\d+)\\)", "\\2", element))
  )

decile_labels <- paste0("Decile ", 3:10)

evo_map <- evo_map |>
  mutate(
    row_label = factor(decile_labels[row + 1], levels = decile_labels),
    col_label = factor(decile_labels[col + 1], levels = decile_labels)
  )

evo_map$is_diag <- evo_map$col_label == evo_map$row_label

evoVCV <- ggplot(evo_map, aes(x = col_label, y = fct_rev(row_label), fill = map)) +
  geom_tile() +
  geom_tile(data = subset(evo_map, is_diag), color = "black", linewidth = 1.5, fill = NA) +
  geom_text(aes(label = round(map, 2)), size = 3, color = "black") +
  scale_fill_gradient(
    low    = "white",
    high   = "#a31e22",
    limits = c(shared_min, shared_max)
  )+
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
# Posterior predictive differences between modern humans and neand --------
hs_preds <- read_rds(paste0(input, "/lc/posteriorPredictive/hsPostPred.rds"))
ne_preds <- read_rds(paste0(input, "lc/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds <- read_rds(paste0(input, "lc/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds <- read_rds(paste0(input, "lc/posteriorPredictive/pantroglodytesPostPred.rds"))
gb_preds <- read_rds(paste0(input, "lc/posteriorPredictive/gorrillaberingeiPostPred.rds"))
gg_preds <- read_rds(paste0(input, "lc/posteriorPredictive/gorillagorillaPostPred.rds"))
pa_preds <- read_rds(paste0(input, "lc/posteriorPredictive/pongoabeliiPostPred.rds"))
ppyg_preds <- read_rds(paste0(input, "lc/posteriorPredictive/pongopygmaeusPostPred.rds"))

trait_labels <- paste0("Decile ", 3:10)

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
              vjust = 1.5, size = 3, color = "grey30") +
    scale_fill_manual(values = color_vec) +
    scale_y_continuous(limits = c(0, 40)) +
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

combined <- (homo + pan) / (gorilla + pongo)
combined
ggsave(paste0(output, "postPred.svg"), plot = combined, width = 14, height = 14)


# Posterior predictive ui2 --------
hs_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/hsPostPred.rds"))
ne_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/neanderthalPostPred.rds"))
pp_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/panpaniscusPostPred.rds"))
pt_preds <- read_rds(paste0(input, "ui2/posteriorPredictive/pantroglodytesPostPred.rds"))

trait_labels <- paste0("Decile ", 3:10)

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
              vjust = 1.5, size = 3, color = "grey30") +
    scale_fill_manual(values = color_vec) +
    scale_y_continuous(limits = c(0, 40)) +
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

combined <- (homo + pan)
ggsave(paste0(output, "postPredUI2.svg"), plot = combined, width = 14, height = 7)

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
hs_preds <- read_rds(paste0(input, "/lc/posteriorPredictive/hsPostPred.rds"))
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

# LC intraspecific means vs. MLE ------------------------------------------

plot_species_posteriors <- function(lc_posterior, lc_mle, species, bins = 500) {
  
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
    scale_x_continuous(limits = c(0,30))+
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



# missing data imputation -------------------------------------------------

lc_dat <- read.csv("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/LCdec3_10.csv")
species <- "Pongo_abelii"
pan_pan_lc <- filter(lc_dat, genus == species)

decile_cols <- c("Decile.3", "Decile.4", "Decile.5", "Decile.6",
                 "Decile.7", "Decile.8", "Decile.9", "Buccal.decile.10..cervical.")

decile_label_map <- c(
  "Decile.3"                    = "Decile 3",
  "Decile.4"                    = "Decile 4",
  "Decile.5"                    = "Decile 5",
  "Decile.6"                    = "Decile 6",
  "Decile.7"                    = "Decile 7",
  "Decile.8"                    = "Decile 8",
  "Decile.9"                    = "Decile 9",
  "Buccal.decile.10..cervical." = "Decile 10"
)

decile_levels <- paste("Decile", 3:10)

df_long <- pan_pan_lc %>%
  mutate(id = row_number()) %>%
  pivot_longer(
    cols      = all_of(decile_cols),
    names_to  = "decile",
    values_to = "value"
  ) %>%
  mutate(
    decile_num = factor(decile_label_map[decile], levels = decile_levels),
    is_missing = is.na(value)
  )

y_to_decile  <- paste("Decile", 3:10)
missing_vars <- grep(paste0("^missing_", species, "_"), names(lc_posterior), value = TRUE)
rows         <- vector("list", length(missing_vars))

for (i in seq_along(missing_vars)) {
  vname   <- missing_vars[i]
  nums    <- as.integer(regmatches(vname, gregexpr("[0-9]+", vname))[[1]])
  obs_idx <- nums[length(nums) - 1]
  dec_idx <- nums[length(nums)]
  
  samps   <- as.numeric(lc_posterior[[vname]])
  map_est <- as.numeric(bayestestR::map_estimate(samps))
  
  rows[[i]] <- data.frame(
    varname    = vname,
    id         = obs_idx + 1L,
    dec_idx    = dec_idx,
    map_est    = map_est,
    post_value = samps
  )
}

posterior_samples_long <- bind_rows(rows) %>%
  mutate(
    decile_num = factor(y_to_decile[dec_idx + 1L], levels = decile_levels),
    id_char    = as.character(id)
  )

map_lookup <- posterior_samples_long %>%
  distinct(id, decile_num, map_est)

missing_ids <- sort(unique(map_lookup$id))
id_colors   <- setNames(
  c("blue"),
  as.character(missing_ids)
)

posterior_samples_long <- posterior_samples_long %>%
  mutate(fill_color = id_colors[as.character(id)])

df_plot <- df_long %>%
  left_join(map_lookup, by = c("id", "decile_num")) %>%
  mutate(
    plot_value = if_else(is_missing, map_est, value),
    line_color = if_else(id %in% missing_ids, as.character(id), "black")
  )

all_levels  <- unique(df_plot$line_color)
color_scale <- ifelse(all_levels == "black", "black", id_colors[all_levels])
names(color_scale) <- all_levels

p1 <- ggplot() +
  geom_violin(
    data     = posterior_samples_long,
    aes(x    = decile_num, y = post_value,
        group = interaction(id, decile_num),
        fill  = id_char),
    alpha    = 0.4,
    color    = NA,
    scale    = "width",
    position = position_identity()
  ) +
  scale_fill_manual(values = id_colors) +
  geom_line(
    data      = df_plot,
    aes(x = decile_num, y = plot_value, group = id, color = line_color),
    linewidth = 0.5
  ) +
  geom_point(
    data  = filter(df_plot, !is_missing),
    aes(x = decile_num, y = plot_value, group = id, color = line_color),
    size  = 0.8
  ) +
  geom_point(
    data  = filter(df_plot, is_missing),
    aes(x = decile_num, y = plot_value, group = id, color = line_color),
    shape = 21,
    size  = 0.8,
    fill  = "white"
  ) +
  scale_color_manual(values = color_scale) +
  labs(x = NULL, y = "Perikymata per millimeter") +
  theme_minimal(base_family = "Georgia") +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )
p1
ggsave(paste0(output, "pongoAbeliiPred.svg"), plot = p1, width = 7, height = 6)

# simulated data posteriors ----------------------------------------------------------

sim_post <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/simulatedData/simulatedDataResults.tsv")
#sim_post <- sim_post[round(0.1 * nrow(sim_post)) : nrow(sim_post), ] #apply burnin
true_missing_val <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/exampleSimulatedData/trueMissingValue.tsv", header = FALSE)[1,1]

#### simulated missing vla trace plot ####
missing_trace <- data.frame(
  idx = sim_post$n,
  post = sim_post$missing_Pan_paniscus_.7.5.)

burnin_cutoff <- 0.1 * max(missing_trace$idx)

p1 <- ggplot(data = missing_trace) +
  geom_point(
    aes(x = idx, 
        y = post,
        color = case_when(
          idx <= burnin_cutoff                ~ "burnin",
          post <= quantile(post, 0.025)       ~ "tail",
          post >= quantile(post, 0.975)       ~ "tail",
          TRUE                               ~ "middle"
        )),
    alpha = 0.05
  ) +
  geom_line(
    aes(x = idx, y = post,
        color = case_when(
          idx <= burnin_cutoff ~ "burnin",
          TRUE                 ~ "middle"
        )),
    alpha = 0.1
  ) +
  geom_hline(
    yintercept = true_missing_val,
    color      = "blue",
    linewidth  = 0.8,
    linetype   = "solid"
  ) +
  scale_color_manual(
    values = c("burnin" = "grey60", "tail" = "red", "middle" = "black"),
    labels = c("burnin" = "Burn-in", "tail" = "Lower/Upper 2.5%", "middle" = "Middle 95%"),
    name   = NULL
  ) +
  labs(
    x = "Cycle",
    y = "Imputed perikymata per millimeter"
  ) +
  theme_minimal(base_family = "Georgia") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )
p1 

ggsave(
  paste0(output, "simulatedMissingDataTrace.pdf"), 
  plot = p1, 
  width = 10, height = 8,
  device = cairo_pdf
)

#### evolutionary vcv ####
#now apply burn in
sim_post <- sim_post[round(0.1 * nrow(sim_post)) : nrow(sim_post), ]

trueEvo <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/exampleSimulatedData/evolutionaryVCV.tsv", header = FALSE)

plot_vcv_posterior <- function(posterior, true_mat, target, bins = 500) {
  p <- 8
  decile_labels <- paste0("Decile ", 3:10)
  
  if (target == "evo") {
    prefix      <- "evo_vcv_."
    plot_title  <- "Evolutionary VCV posterior distributions"
    col_builder <- function(i, j) paste0("evo_vcv_.", i, ".", j, ".")
  } else {
    prefix      <- paste0(target, "_vcv_.")
    plot_title  <- bquote("Intraspecific VCV posterior distributions —" ~
                            italic(.(gsub("_", " ", target))))
    col_builder <- function(i, j) paste0(target, "_vcv_.", i, ".", j, ".")
  }
  
  vcv_cols <- as.vector(outer(0:(p-1), 0:(p-1), col_builder))
  missing_cols <- setdiff(vcv_cols, colnames(posterior))
  if (length(missing_cols) > 0)
    stop(sprintf("Missing columns for '%s': %s",
                 target, paste(missing_cols, collapse = ", ")))
  
  vcv_long <- posterior[, vcv_cols, drop = FALSE] %>%
    pivot_longer(cols      = everything(),
                 names_to  = "element",
                 values_to = "value") %>%
    mutate(
      row     = as.integer(sub(paste0(gsub("\\.", "\\\\.", prefix), "(\\d+)\\.(\\d+)\\."), "\\1", element)),
      col     = as.integer(sub(paste0(gsub("\\.", "\\\\.", prefix), "(\\d+)\\.(\\d+)\\."), "\\2", element)),
      row_lbl = factor(decile_labels[row + 1], levels = decile_labels),
      col_lbl = factor(decile_labels[col + 1], levels = decile_labels)
    )
  
  vcv_quantiles <- vcv_long %>%
    group_by(element) %>%
    summarise(
      q_lo   = quantile(value, 0.025),
      q_hi   = quantile(value, 0.975),
      trunc_lo = quantile(value, 0.005),   # <-- new: truncation bounds
      trunc_hi = quantile(value, 0.995),   # <-- new
      .groups = "drop"
    )
  
  true_long <- data.frame(
    row      = rep(0:(p-1), each = p),
    col      = rep(0:(p-1), times = p),
    true_val = as.vector(as.matrix(true_mat))
  ) %>%
    mutate(
      row_lbl = factor(decile_labels[row + 1], levels = decile_labels),
      col_lbl = factor(decile_labels[col + 1], levels = decile_labels)
    )
  
  vcv_long <- vcv_long %>%
    left_join(vcv_quantiles, by = "element") %>%
    mutate(fill_color = if_else(value < q_lo | value > q_hi, "tail", "middle")) %>%
    filter(value >= trunc_lo & value <= trunc_hi)               # <-- truncate here
  
  ggplot(vcv_long, aes(x = value, fill = fill_color)) +
    geom_histogram(bins = bins, color = NA) +
    geom_vline(
      data      = true_long,
      aes(xintercept = true_val),
      color     = "blue",
      linewidth = 0.6,
      linetype  = "solid"
    ) +
    scale_fill_manual(
      values = c("tail" = "red", "middle" = "black"),
      labels = c("tail" = "Lower/Upper 2.5%", "middle" = "Middle 95%"),
      name   = NULL
    ) +
    guides(
      fill = guide_legend(override.aes = list(alpha = 1))
    ) +
    facet_grid(row_lbl ~ col_lbl, scales = "free") +
    labs(
      x     = "Posterior VCV value",
      y     = "Posterior sample count",
      title = plot_title
    ) +
    theme_minimal(base_family = "Georgia") +
    theme(
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA),  # <-- outline
      legend.position  = "none",
      strip.text       = element_text(size = 7),
      axis.text        = element_text(size = 6),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
}

p1 <- plot_vcv_posterior(sim_post, trueEvo, target = "evo")
p1
ggsave(paste0(output, "simulatedEvoVCV.svg"), plot = p1, width = 6, height = 6)
ggsave(
  paste0(output, "simulatedEvoVCV.pdf"), 
  plot = p1, 
  width = 6, height = 6,
  device = cairo_pdf
)


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
