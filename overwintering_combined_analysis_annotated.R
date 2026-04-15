# ============================================================
# Overwintering queens: combined annotated analysis script
# ============================================================
#
# Purpose:
# Recreate the analyses from the prior two analysis runs using
# only the uploaded files:
#   1. overwintering_data_clean.csv
#   2. beeDist.csv
#
# This script includes:
# - mass by season
# - pathogen prevalence by season
# - positive-only pathogen load comparisons
# - positive-only load models with bee mass as a covariate
# - exploratory weight-by-infection-status comparisons
# - fall-only hibernaculum structure summaries and figures
# - multivariate PCA + LDA using mass and pathogen loads
# - fall geometry/network figure based on Euclidean distance
#
# Notes:
# - The pathogen file is already merged and cleaned.
# - The spatial file is fall-only.
# - No GitHub inputs are used anywhere in this script.
# - Place this script in the same working directory as the two CSV files.
#
# Optional package install block:
# install.packages(c("readr", "dplyr", "tidyr", "ggplot2", "broom", "knitr", "MASS", "tibble", "scales"))
# ============================================================

# -----------------------------
# Load packages
# -----------------------------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(knitr)
library(MASS)
library(tibble)
library(scales)

theme_set(theme_bw(base_size = 12))
options(dplyr.summarise.inform = FALSE)

# -----------------------------
# Import and prepare data
# -----------------------------
# The pathogen data are already in long format:
# one row per queen-pathogen combination.
pathogen_df <- read_csv("overwintering_data_clean.csv", show_col_types = FALSE) %>%
  rename(
    sample_id = ID,
    season = Season,
    pathogen = pathogen,
    load = load,
    weight_g = `Bee_Mass_Est_g.`
  ) %>%
  mutate(
    season = factor(season, levels = c("Fall", "Spring")),
    positive = replace_na(load, 0) > 0,
    log10_load = if_else(load > 0, log10(load), NA_real_),
    log10_load_plus1 = log10(replace_na(load, 0) + 1)
  ) %>%
  select(sample_id, season, pathogen, load, weight_g, positive, log10_load, log10_load_plus1, time)

# One row per queen for season-level mass analyses.
weight_df <- pathogen_df %>%
  distinct(sample_id, season, weight_g)

# Fall-only spatial dataset.
# The weight column in this file is not used below because the merged
# pathogen file already contains the mass values used in the models.
fall_holes <- read_csv("beeDist.csv", show_col_types = FALSE) %>%
  rename(
    hole = Hole,
    nest_id = ID,
    lab_id = LabID,
    surface_distance_cm = distance_from_hole_on_surface_cm,
    nearest_neighbor_cm = nearest_neighbor
  ) %>%
  mutate(
    hole = factor(hole),
    euclidean_cm = sqrt(depth_cm^2 + surface_distance_cm^2)
  )

# Helper for exact binomial confidence intervals on prevalence.
exact_ci <- function(x, n) {
  bt <- binom.test(x, n)
  tibble(ci_low = bt$conf.int[1], ci_high = bt$conf.int[2])
}

# Quick console check so the user can confirm file dimensions.
cat("\n================ DATA OVERVIEW ================\n")
cat("Merged pathogen rows:", nrow(pathogen_df), "\n")
cat("Unique queens:", n_distinct(pathogen_df$sample_id), "\n")
cat("Unique pathogens:", n_distinct(pathogen_df$pathogen), "\n")
cat("Fall spatial rows:", nrow(fall_holes), "\n")
cat("Occupied holes:", n_distinct(fall_holes$hole), "\n")

# -----------------------------
# 1. Bee mass by season
# -----------------------------
cat("\n================ 1. MASS BY SEASON ================\n")

weight_summary <- weight_df %>%
  group_by(season) %>%
  summarise(
    n = sum(!is.na(weight_g)),
    mean_g = mean(weight_g, na.rm = TRUE),
    sd_g = sd(weight_g, na.rm = TRUE),
    median_g = median(weight_g, na.rm = TRUE),
    min_g = min(weight_g, na.rm = TRUE),
    max_g = max(weight_g, na.rm = TRUE)
  )

weight_welch <- t.test(weight_g ~ season, data = weight_df)
weight_wilcox <- wilcox.test(weight_g ~ season, data = weight_df, exact = FALSE)

print(kable(weight_summary, digits = 3))
print(kable(
  tibble(
    comparison = "Fall vs Spring mass",
    test = c("Welch t-test", "Wilcoxon rank-sum"),
    statistic = c(unname(weight_welch$statistic), unname(weight_wilcox$statistic)),
    p_value = c(weight_welch$p.value, weight_wilcox$p.value)
  ),
  digits = 3
))

p_weight <- ggplot(weight_df, aes(season, weight_g, fill = season)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Fall" = "#D4A017", "Spring" = "#4C9A2A")) +
  labs(x = NULL, y = "Bee mass (g)", title = "Bee mass by season") +
  theme(legend.position = "none")
print(p_weight)

# -----------------------------
# 2. Pathogen prevalence by season
# -----------------------------
cat("\n================ 2. PREVALENCE BY SEASON ================\n")

prevalence_summary <- pathogen_df %>%
  group_by(pathogen, season) %>%
  summarise(
    positive_n = sum(positive),
    total_n = n(),
    prevalence = positive_n / total_n,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(ci = list(exact_ci(positive_n, total_n))) %>%
  unnest_wider(ci) %>%
  ungroup()

prevalence_tests <- pathogen_df %>%
  group_by(pathogen) %>%
  group_modify(~{
    tab <- table(.x$season, .x$positive)
    ft <- fisher.test(tab)
    tibble(
      fall_positive = sum(.x$positive[.x$season == "Fall"]),
      fall_n = sum(.x$season == "Fall"),
      spring_positive = sum(.x$positive[.x$season == "Spring"]),
      spring_n = sum(.x$season == "Spring"),
      odds_ratio = unname(ft$estimate),
      fisher_p = ft$p.value
    )
  }) %>%
  ungroup()

print(kable(prevalence_summary, digits = 3))
print(kable(prevalence_tests, digits = 3))

p_prevalence <- ggplot(prevalence_summary, aes(pathogen, prevalence, color = season)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = position_dodge(width = 0.3),
    width = 0.12,
    linewidth = 0.5
  ) +
  scale_color_manual(values = c("Fall" = "#D4A017", "Spring" = "#4C9A2A")) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(
    x = NULL,
    y = "Prevalence",
    title = "Pathogen prevalence by season (95% exact CI)"
  ) +
  theme(legend.title = element_blank())
print(p_prevalence)

# -----------------------------
# 3. Positive-only pathogen loads
# -----------------------------
cat("\n================ 3. POSITIVE-ONLY LOADS ================\n")

positive_load_summary <- pathogen_df %>%
  filter(positive) %>%
  group_by(pathogen, season) %>%
  summarise(
    n_positive = n(),
    mean_log10 = mean(log10_load, na.rm = TRUE),
    median_log10 = median(log10_load, na.rm = TRUE),
    .groups = "drop"
  )

load_tests <- pathogen_df %>%
  filter(positive) %>%
  group_by(pathogen) %>%
  group_modify(~{
    fall_vals <- .x %>% filter(season == "Fall") %>% pull(log10_load) %>% na.omit()
    spring_vals <- .x %>% filter(season == "Spring") %>% pull(log10_load) %>% na.omit()

    welch_t <- if (length(fall_vals) > 1 && length(spring_vals) > 1) {
      unname(t.test(fall_vals, spring_vals)$statistic)
    } else {
      NA_real_
    }

    welch_p <- if (length(fall_vals) > 1 && length(spring_vals) > 1) {
      t.test(fall_vals, spring_vals)$p.value
    } else {
      NA_real_
    }

    wilcoxon_u <- if (length(fall_vals) > 0 && length(spring_vals) > 0) {
      unname(wilcox.test(fall_vals, spring_vals, exact = FALSE)$statistic)
    } else {
      NA_real_
    }

    wilcoxon_p <- if (length(fall_vals) > 0 && length(spring_vals) > 0) {
      wilcox.test(fall_vals, spring_vals, exact = FALSE)$p.value
    } else {
      NA_real_
    }

    tibble(
      n_fall_positive = length(fall_vals),
      n_spring_positive = length(spring_vals),
      welch_t = welch_t,
      welch_p = welch_p,
      wilcoxon_u = wilcoxon_u,
      wilcoxon_p = wilcoxon_p
    )
  }) %>%
  ungroup()

print(kable(positive_load_summary, digits = 3))
print(kable(load_tests, digits = 3))

p_load <- ggplot(pathogen_df %>% filter(positive), aes(season, log10_load, fill = season)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 1.8, alpha = 0.75) +
  scale_fill_manual(values = c("Fall" = "#D4A017", "Spring" = "#4C9A2A")) +
  facet_wrap(~ pathogen, scales = "free_y") +
  labs(x = NULL, y = "log10(load)", title = "Positive-only pathogen loads by season") +
  theme(legend.position = "none")
print(p_load)

# -----------------------------
# 4. Positive-only load models with bee mass
# -----------------------------
cat("\n================ 4. LOAD MODELS WITH MASS COVARIATE ================\n")

positive_load_models <- pathogen_df %>%
  filter(pathogen %in% c("BQCV", "Nosema"), positive) %>%
  group_by(pathogen) %>%
  group_modify(~{
    fit <- lm(log10_load ~ season + weight_g, data = .x)
    broom::tidy(fit) %>%
      mutate(
        n = nrow(.x),
        r_squared = summary(fit)$r.squared
      )
  }) %>%
  ungroup()

print(kable(positive_load_models, digits = 3))

p_weight_load <- ggplot(
  pathogen_df %>% filter(pathogen %in% c("BQCV", "Nosema"), positive),
  aes(weight_g, log10_load, color = season)
) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_color_manual(values = c("Fall" = "#D4A017", "Spring" = "#4C9A2A")) +
  facet_wrap(~ pathogen, scales = "free_y") +
  labs(
    x = "Bee mass (g)",
    y = "log10(load)",
    title = "Positive-only load versus bee mass"
  ) +
  theme(legend.title = element_blank())
print(p_weight_load)

# -----------------------------
# 5. Exploratory weight by infection status
# -----------------------------
cat("\n================ 5. WEIGHT BY INFECTION STATUS ================\n")

weight_by_infection <- pathogen_df %>%
  group_by(pathogen) %>%
  group_modify(~{
    pos_w <- .x %>% filter(positive) %>% pull(weight_g) %>% na.omit()
    neg_w <- .x %>% filter(!positive) %>% pull(weight_g) %>% na.omit()

    welch_t <- if (length(pos_w) > 1 && length(neg_w) > 1) {
      unname(t.test(pos_w, neg_w)$statistic)
    } else {
      NA_real_
    }

    welch_p <- if (length(pos_w) > 1 && length(neg_w) > 1) {
      t.test(pos_w, neg_w)$p.value
    } else {
      NA_real_
    }

    wilcox_u <- if (length(pos_w) > 1 && length(neg_w) > 1) {
      unname(wilcox.test(pos_w, neg_w, exact = FALSE)$statistic)
    } else {
      NA_real_
    }

    wilcox_p <- if (length(pos_w) > 1 && length(neg_w) > 1) {
      wilcox.test(pos_w, neg_w, exact = FALSE)$p.value
    } else {
      NA_real_
    }

    tibble(
      n_positive = length(pos_w),
      n_negative = length(neg_w),
      mean_weight_positive_g = mean(pos_w),
      mean_weight_negative_g = mean(neg_w),
      welch_t = welch_t,
      welch_p = welch_p,
      wilcoxon_u = wilcox_u,
      wilcoxon_p = wilcox_p
    )
  }) %>%
  ungroup()

print(kable(weight_by_infection, digits = 3))

# -----------------------------
# 6. Fall-only hibernaculum structure
# -----------------------------
cat("\n================ 6. FALL HIBERNACULUM STRUCTURE ================\n")

hole_occupancy <- fall_holes %>%
  count(hole, name = "queens_in_hole")

fall_spatial_summary <- tibble(
  metric = c(
    "Queens in spatial dataset",
    "Occupied holes",
    "Mean queens per occupied hole",
    "Maximum queens in one hole",
    "Mean depth (cm)",
    "Median depth (cm)",
    "Mean surface distance (cm)",
    "Nearest-neighbor observations"
  ),
  value = c(
    nrow(fall_holes),
    n_distinct(fall_holes$hole),
    mean(hole_occupancy$queens_in_hole),
    max(hole_occupancy$queens_in_hole),
    mean(fall_holes$depth_cm, na.rm = TRUE),
    median(fall_holes$depth_cm, na.rm = TRUE),
    mean(fall_holes$surface_distance_cm, na.rm = TRUE),
    sum(!is.na(fall_holes$nearest_neighbor_cm))
  )
)

fall_spatial_correlations <- tibble(
  comparison = c(
    "Depth vs surface distance",
    "Depth vs nearest-neighbor distance",
    "Surface distance vs nearest-neighbor distance"
  ),
  rho = c(
    cor(fall_holes$depth_cm, fall_holes$surface_distance_cm, method = "spearman", use = "complete.obs"),
    cor(fall_holes$depth_cm, fall_holes$nearest_neighbor_cm, method = "spearman", use = "complete.obs"),
    cor(fall_holes$surface_distance_cm, fall_holes$nearest_neighbor_cm, method = "spearman", use = "complete.obs")
  ),
  p_value = c(
    cor.test(fall_holes$depth_cm, fall_holes$surface_distance_cm, method = "spearman", exact = FALSE)$p.value,
    cor.test(fall_holes$depth_cm, fall_holes$nearest_neighbor_cm, method = "spearman", exact = FALSE)$p.value,
    cor.test(fall_holes$surface_distance_cm, fall_holes$nearest_neighbor_cm, method = "spearman", exact = FALSE)$p.value
  )
)

print(kable(fall_spatial_summary, digits = 3))
print(kable(hole_occupancy, digits = 3))
print(kable(fall_spatial_correlations, digits = 3))

p_depth_hist <- ggplot(fall_holes, aes(depth_cm)) +
  geom_histogram(bins = 7, fill = "#4C78A8", color = "white") +
  labs(
    x = "Depth under soil (cm)",
    y = "Frequency",
    title = "Fall hibernaculum depth distribution"
  )
print(p_depth_hist)

p_hole_occupancy <- ggplot(hole_occupancy, aes(hole, queens_in_hole)) +
  geom_col(fill = "#4C78A8") +
  labs(
    x = "Hole",
    y = "Queens recovered",
    title = "Queens recovered per occupied fall hole"
  )
print(p_hole_occupancy)

p_depth_surface <- ggplot(fall_holes, aes(surface_distance_cm, depth_cm)) +
  geom_point(size = 2.3, color = "#4C78A8") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  labs(
    x = "Distance from surface hole (cm)",
    y = "Depth (cm)",
    title = "Depth vs surface distance"
  )
print(p_depth_surface)

p_depth_neighbor <- ggplot(
  fall_holes %>% filter(!is.na(nearest_neighbor_cm)),
  aes(nearest_neighbor_cm, depth_cm)
) +
  geom_point(size = 2.3, color = "#4C78A8") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  labs(
    x = "Nearest-neighbor distance (cm)",
    y = "Depth (cm)",
    title = "Depth vs nearest-neighbor distance"
  )
print(p_depth_neighbor)

# -----------------------------
# 7. Multivariate season prediction
# -----------------------------
cat("\n================ 7. MULTIVARIATE PCA + LDA ================\n")

# Convert the long pathogen table to one row per queen.
multivar_df <- pathogen_df %>%
  mutate(load = replace_na(load, 0)) %>%
  select(sample_id, season, pathogen, load, weight_g) %>%
  pivot_wider(names_from = pathogen, values_from = load) %>%
  mutate(
    log10_BQCV = log10(BQCV + 1),
    log10_DWV = log10(DWV + 1),
    log10_Nosema = log10(Nosema + 1)
  ) %>%
  select(sample_id, season, weight_g, log10_BQCV, log10_DWV, log10_Nosema) %>%
  drop_na()

multivariate_summary <- multivar_df %>%
  group_by(season) %>%
  summarise(
    n = n(),
    mean_weight_g = mean(weight_g),
    mean_log10_BQCV = mean(log10_BQCV),
    mean_log10_DWV = mean(log10_DWV),
    mean_log10_Nosema = mean(log10_Nosema)
  )
print(kable(multivariate_summary, digits = 3))

# PCA on centered/scaled predictors.
pca_fit <- prcomp(
  multivar_df %>% select(weight_g, log10_BQCV, log10_DWV, log10_Nosema),
  center = TRUE,
  scale. = TRUE
)

scores <- as_tibble(pca_fit$x[, 1:2]) %>%
  bind_cols(multivar_df %>% select(sample_id, season))

loadings <- as.data.frame(pca_fit$rotation[, 1:2]) %>%
  rownames_to_column("variable")

var_explained <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
arrow_scale <- 1.8

loadings_plot <- loadings %>%
  mutate(
    PC1 = PC1 * arrow_scale,
    PC2 = PC2 * arrow_scale
  )

season_centroids <- scores %>%
  group_by(season) %>%
  summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2))

print(kable(loadings, digits = 3))
print(kable(season_centroids, digits = 3))

p_pca <- ggplot(scores, aes(PC1, PC2, color = season)) +
  stat_ellipse(level = 0.68, linewidth = 0.7, alpha = 0.6, show.legend = FALSE) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_segment(
    data = loadings_plot,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.18, "cm")),
    linewidth = 0.7,
    color = "gray30"
  ) +
  geom_text(
    data = loadings_plot,
    aes(x = PC1, y = PC2, label = variable),
    inherit.aes = FALSE,
    size = 3.6,
    color = "gray20",
    nudge_x = 0.05,
    nudge_y = 0.05
  ) +
  scale_color_manual(values = c("Fall" = "#D4A017", "Spring" = "#4C9A2A")) +
  labs(
    title = "PCA of queen mass and pathogen loads",
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "% variance)"),
    color = NULL
  )
print(p_pca)

# LDA provides an exploratory classification summary.
lda_fit <- MASS::lda(
  season ~ weight_g + log10_BQCV + log10_DWV + log10_Nosema,
  data = multivar_df
)

lda_cv <- MASS::lda(
  season ~ weight_g + log10_BQCV + log10_DWV + log10_Nosema,
  data = multivar_df,
  CV = TRUE
)

conf_mat <- table(
  Observed = multivar_df$season,
  Predicted = lda_cv$class
)

accuracy <- mean(lda_cv$class == multivar_df$season)

lda_coefficients <- broom::tidy(lda_fit) %>%
  select(term, LD1)

print(kable(lda_coefficients, digits = 3))
print(kable(as.data.frame.matrix(conf_mat) %>% rownames_to_column("Observed"), digits = 3))
print(kable(tibble(metric = "Leave-one-out cross-validated accuracy", value = accuracy), digits = 3))

# -----------------------------
# 8. Fall geometry / network figure
# -----------------------------
cat("\n================ 8. FALL GEOMETRY NETWORK ================\n")

fall_distance_summary <- fall_holes %>%
  summarise(
    n_queens = n(),
    n_holes = n_distinct(hole),
    mean_depth_cm = mean(depth_cm, na.rm = TRUE),
    mean_surface_distance_cm = mean(surface_distance_cm, na.rm = TRUE),
    mean_euclidean_cm = mean(euclidean_cm, na.rm = TRUE),
    median_euclidean_cm = median(euclidean_cm, na.rm = TRUE),
    max_euclidean_cm = max(euclidean_cm, na.rm = TRUE)
  )

by_hole_distance <- fall_holes %>%
  group_by(hole) %>%
  summarise(
    queens_in_hole = n(),
    mean_depth_cm = mean(depth_cm),
    mean_surface_distance_cm = mean(surface_distance_cm),
    mean_euclidean_cm = mean(euclidean_cm)
  )

print(kable(fall_distance_summary, digits = 3))
print(kable(by_hole_distance, digits = 3))

hole_nodes <- fall_holes %>%
  distinct(hole) %>%
  mutate(x = 0, y = 0)

p_geometry_network <- ggplot(fall_holes) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = surface_distance_cm,
      yend = -depth_cm,
      linewidth = euclidean_cm
    ),
    color = "gray35",
    alpha = 0.8
  ) +
  geom_point(
    data = hole_nodes,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    size = 3.6,
    color = "black"
  ) +
  geom_point(
    aes(x = surface_distance_cm, y = -depth_cm),
    size = 2.8,
    color = "#2C7FB8"
  ) +
  geom_text(
    aes(x = surface_distance_cm, y = -depth_cm, label = lab_id),
    nudge_y = -0.45,
    size = 2.8,
    check_overlap = TRUE
  ) +
  scale_linewidth_continuous(range = c(0.4, 2.2)) +
  coord_equal() +
  facet_wrap(~ hole, scales = "fixed") +
  labs(
    title = "Fall queen subsurface geometry by focal hole",
    subtitle = "Each edge runs from the focal surface hole (black node) to an excavated queen; edge width encodes Euclidean distance.",
    x = "Horizontal distance from focal hole (cm)",
    y = "Depth below surface (cm)",
    linewidth = "Euclidean\ndistance (cm)"
  )
print(p_geometry_network)

# -----------------------------
# Optional save lines
# -----------------------------
# Uncomment any of these if you want written figure files.
# ggsave("fig_weight_by_season.png", p_weight, width = 6.8, height = 4.5, dpi = 300)
# ggsave("fig_prevalence_by_season.png", p_prevalence, width = 7.2, height = 4.8, dpi = 300)
# ggsave("fig_positive_only_loads.png", p_load, width = 10.5, height = 4.6, dpi = 300)
# ggsave("fig_weight_covariate_loads.png", p_weight_load, width = 9.6, height = 4.4, dpi = 300)
# ggsave("fig_fall_depth_histogram.png", p_depth_hist, width = 6.8, height = 4.4, dpi = 300)
# ggsave("fig_fall_hole_occupancy.png", p_hole_occupancy, width = 6.8, height = 4.4, dpi = 300)
# ggsave("fig_depth_vs_surface_distance.png", p_depth_surface, width = 6.8, height = 4.4, dpi = 300)
# ggsave("fig_depth_vs_neighbor_distance.png", p_depth_neighbor, width = 6.8, height = 4.4, dpi = 300)
# ggsave("fig_pca_multivariate.png", p_pca, width = 7.8, height = 6.4, dpi = 300)
# ggsave("fig_fall_geometry_network.png", p_geometry_network, width = 11, height = 7.6, dpi = 300)

# -----------------------------
# Session info for reproducibility
# -----------------------------
cat("\n================ SESSION INFO ================\n")
print(sessionInfo())
