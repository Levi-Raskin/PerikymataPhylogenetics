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


# Evo VCV combined --------------------------------------------------------

map_estimate <- function(x) {
  as.numeric(bayestestR::map_estimate(x))
}
#evo VCV MAP
evo_vcv_cols <- paste0("evo_vcv_(", rep(0:7, each = 8), ",", rep(0:7, times = 8), ")")

shared_max <- max(
  dplyr::select(lc_posterior, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist(),
  dplyr::select(lc_posterior_no_hominin, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist(),
  dplyr::select(ui2_posterior, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist()
)

shared_min <- min(
  dplyr::select(lc_posterior, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist(),
  dplyr::select(lc_posterior_no_hominin, all_of(evo_vcv_cols)) |> 
    dplyr::summarise(across(everything(), map_estimate)) |> unlist(),
  dplyr::select(ui2_posterior, all_of(evo_vcv_cols)) |> 
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

lcHominins <- treeplot +  evoVCV
lcHominins

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

lcNoHominins <- treeplot + evoVCV
lcNoHominins

# UI2 evo VCV ---------------------------------------------------------------

evo_vcv_cols <- paste0("evo_vcv_(", rep(0:7, each = 8), ",", rep(0:7, times = 8), ")")

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

plottree <- ape::read.tree(file = "/Users/levir/Documents/GitHub/PerikymataPhylogenetics/data/tree.txt")
plottree <- keep.tip(plottree, c("Neanderthal", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"))
plottree$tip.label <- gsub("_", " ", plottree$tip.label)
treeplot <- ggtree(plottree) + 
  geom_tiplab(aes(fontface = ifelse(label %in% c("Modern humans", "Neanderthal"), 2, 4)), family = "Georgia") +
  hexpand(0.55)
treeplot <- ggtree::rotate(treeplot, 5)
treeplot

ui2 <- treeplot + evoVCV
ui2

combined <- lcHominins / lcNoHominins / ui2
combined
ggsave(paste0(output, "treeMAPCombined.svg"), plot = combined, width = 14, height = 12)

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
      legend.position = "right",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
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
    c(sapply(trait_labels, function(d) TeX(paste0(d, " — $\\bar{C}$"))),
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
ggsave(paste0(output, "postPredCombined.svg"), plot = combined, width = 14, height = 21)


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
    rename(any_of(c("Decile.10" = "Buccal.decile.10..cervical.")))
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
    rename(any_of(c("Decile.10" = "Buccal.decile.10..cervical.")))
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

#LC
lc_mean_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_posterior_predictive_check_means.RDS")
lc_sd_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/lc/lc_posterior_predictive_check_sd.RDS")

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
ui2_mean_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_posterior_predictive_check_means.RDS")
ui2_sd_list <- readRDS("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/withGibbs_v2/ui2/ui2_posterior_predictive_check_sd.RDS")

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

total <- matrix(data = NA, nrow = 0, ncol = 8)
for (i in species) {
  p <- calc_ppc_pvalues(
    ppc_means    = ui2_mean_list[[i]],
    ppc_sds      = ui2_sd_list[[i]],
    observed_data = ui2_data,
    species      = i
  )
  total <- rbind(total, p)
}
