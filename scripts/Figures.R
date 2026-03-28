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

lc_posterior <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8.tsv")
lc_posterior <- lc_posterior[round(0.1 * nrow(lc_posterior)) : nrow(lc_posterior), ] #apply burnin

lc_posterior_no_hominin <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8_no_hominin.tsv")
lc_posterior_no_hominin  <- lc_posterior_no_hominin[round(0.1 * nrow(lc_posterior_no_hominin)) : nrow(lc_posterior_no_hominin), ] #apply burnin

lc_posterior_species_means <- read.delim("/Users/levir/Documents/GitHub/PerikymataPhylogenetics/results/lc/lc_dec3_8_species_means.tsv")
lc_posterior_species_means  <- lc_posterior_species_means[round(0.1 * nrow(lc_posterior_species_means)) : nrow(lc_posterior_species_means), ] #apply burnin

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

#evo VCV MAP
evo_vcv_cols <- paste0("evo_vcv_.", rep(0:7, each = 8), ".", rep(0:7, times = 8), ".")

evo_map <- dplyr::select(lc_posterior, all_of(evo_vcv_cols)) |>
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

n_samples <- nrow(lc_posterior)
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
# hs_preds <- get_posterior_predictive(lc_posterior, "Homo_sapiens", n_traits, n_samples)
# ne_preds <- get_posterior_predictive(lc_posterior, "Neanderthal", n_traits, n_samples)
# pp_preds <- get_posterior_predictive(lc_posterior, "Pan_paniscus", n_traits, n_samples)
# pt_preds <- get_posterior_predictive(lc_posterior, "Pan_troglodytes", n_traits, n_samples)
# gb_preds <- get_posterior_predictive(lc_posterior, "Gorilla_beringei", n_traits, n_samples)
# gg_preds <- get_posterior_predictive(lc_posterior, "Gorilla_gorilla", n_traits, n_samples)
# pa_preds <- get_posterior_predictive(lc_posterior, "Pongo_abelii", n_traits, n_samples)
# ppyg_preds <- get_posterior_predictive(lc_posterior, "Pongo_pygmaeus", n_traits, n_samples)
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
  
  if (!species %in% lc_mle$genus) {
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
    mutate(decile = factor(decile_map[mean_col], levels = decile_levels))
  
  sp_quantiles <- sp_long %>%
    group_by(decile) %>%
    summarise(q_lo = quantile(value, 0.025),
              q_hi = quantile(value, 0.975),
              .groups = "drop")
  
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
    scale_x_continuous(limits = c(0,50))+
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
    paste0(output, "/meanPosteriorHists/",i,".svg"), 
    plot = p, 
    width = 10, height = 8
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
