## Plots of the significant organelle x organelle associations from 04_organelle_associations.R
## 1) association heatmap   2) odds ratios + frequency shift   3) contingency tables
## Output: plots/assoc_organelles_*

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- NA_character_   # NA = read the main sheet from Google Sheets; or a tsv snapshot
default_outdir     <- "C:/projects/tem_classification"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/tem_classification/analysis"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "assoc_common.R"))
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)

cor_data <- build_cor_data(load_dino_data(opt))
to_be_plotted <- read_assoc_results("assoc_organelles", opt$table_dir, "04_organelle_associations.R") %>%
  filter(p_adj < 0.05)
if (nrow(to_be_plotted) == 0) {
  msg <- "No significant organelle x organelle pairs, nothing to plot"
  if (interactive()) stop(msg) else { message(msg); quit(save = "no") }
}

## ---- 1) association heatmap: every significant pair ----
plot_pairs <- to_be_plotted %>% filter(!is.na(estimate))
feats <- unique(c(plot_pairs$var1, plot_pairs$var2))
m <- matrix(0, length(feats), length(feats), dimnames = list(feats, feats))
m[cbind(plot_pairs$var1, plot_pairs$var2)] <- plot_pairs$estimate
m <- m + t(m)
feat_clustered <- if (length(feats) > 2) feats[hclust(dist(m))$order] else feats

## manual order: chloroplast related, common organelles, specific ones (clustered), undescribed at the end
undesc_feats <- str_sort(str_subset(feats, "^undescribed_"), numeric = TRUE)
specific_feats <- setdiff(feat_clustered, c(chloro_feats, common_feats, undesc_feats))
feat_order <- c(intersect(feat_clustered, chloro_feats),
                intersect(common_feats, feats),
                specific_feats,
                undesc_feats)

## lower triangle only: orient every pair so the x feature comes before the y feature
heatmap_data <- plot_pairs %>%
  mutate(i1 = match(var1, feat_order), i2 = match(var2, feat_order),
         x = feat_order[pmin(i1, i2)], y = feat_order[pmax(i1, i2)],
         x_group = feature_group(x), y_group = feature_group(y),
         x = factor(x, feat_order[-length(feat_order)]),
         y = factor(y, rev(feat_order[-1])))

## number of cells with the feature, drawn outside the outer panels:
## x counts above the top facet row, y counts right of the last column
feat_counts <- map_int(cor_data, ~ if (is_binary(.x)) sum(.x == 1, na.rm = TRUE) else sum(!is.na(.x)))
top_counts <- distinct(heatmap_data, x, x_group) %>%
  mutate(y_group = factor(levels(droplevels(heatmap_data$y_group))[1], levels(heatmap_data$y_group)),
         n = feat_counts[as.character(x)])
right_counts <- distinct(heatmap_data, y, y_group) %>%
  mutate(x_group = factor(last(levels(droplevels(heatmap_data$x_group))), levels(heatmap_data$x_group)),
         n = feat_counts[as.character(y)])

## light background only for the group blocks that contain pairs, empty blocks stay invisible
heatmap_blocks <- distinct(heatmap_data, x_group, y_group)

p_heatmap <- ggplot(heatmap_data, aes(x, y, fill = estimate)) +
  geom_rect(data = heatmap_blocks, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = "grey95", inherit.aes = FALSE) +
  geom_tile(colour = "white") +
  geom_text(aes(label = stars(p_adj)), size = 3, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "phi") +
  geom_text(data = top_counts, aes(x = x, y = Inf, label = n), inherit.aes = FALSE,
            angle = 90, hjust = -0.2, size = 2.5, colour = "grey40") +
  geom_text(data = right_counts, aes(x = Inf, y = y, label = n), inherit.aes = FALSE,
            hjust = -0.3, size = 2.5, colour = "grey40") +
  facet_grid(y_group ~ x_group, scales = "free", space = "free", switch = "both") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = "Significant organelle x organelle associations",
       subtitle = "* p_adj < 0.05, ** < 0.01, *** < 0.001 | top / right: number of cells with the feature") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        plot.margin = margin(10, 10, 5, 5), plot.subtitle = element_text(margin = margin(b = 30)),
        legend.box.spacing = unit(30, "pt"),
        panel.spacing = unit(4, "pt"), strip.placement = "outside",
        strip.text = element_text(face = "bold"), strip.text.y.left = element_text(angle = 0, hjust = 1))

save_plot(p_heatmap, "assoc_organelles_heatmap", 10, 9, opt$plot_dir)

## ---- 2) odds ratios + how var1 shifts the frequency of var2 ----
or_cap <- 6   # log2 OR axis limit, OR = 0 / Inf are drawn at the cap
or_pairs <- to_be_plotted %>%
  filter(test == "fisher 2x2") %>%
  separate_wider_delim(var1_present, "/", names = c("a", "b"), cols_remove = FALSE) %>%
  separate_wider_delim(var1_absent, "/", names = c("c", "d"), cols_remove = FALSE) %>%
  mutate(across(c(a, b, c, d), as.numeric),
         prev_var1_present = a / (a + b),
         prev_var1_absent  = c / (c + d),
         pair = fct_reorder(paste(var1, "->", var2), log2(odds_ratio)),
         direction = if_else(odds_ratio > 1, "co-occur", "exclude"),
         point = if_else(abs(log2(odds_ratio)) >= or_cap, "capped (OR 0 / Inf)", "estimate"),
         across(c(odds_ratio, or_ci_low, or_ci_high), ~ pmin(pmax(log2(.x), -or_cap), or_cap), .names = "l2_{.col}"))

if (nrow(or_pairs) > 0) {
  p_forest <- ggplot(or_pairs, aes(l2_odds_ratio, pair, colour = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = l2_or_ci_low, xmax = l2_or_ci_high), width = 0.3, orientation = "y") +
    geom_point(aes(shape = point), size = 2.5) +
    geom_text(aes(x = or_cap + 0.5, label = paste(var1_present, "vs", var1_absent)),
              hjust = 0, colour = "grey30", size = 3) +
    scale_colour_manual(values = c("co-occur" = "#b2182b", "exclude" = "#2166ac"), name = NULL) +
    scale_shape_manual(values = c("estimate" = 16, "capped (OR 0 / Inf)" = 4), name = NULL) +
    coord_cartesian(xlim = c(-or_cap, or_cap + 3), clip = "off") +
    labs(x = "log2 odds ratio (95% CI)", y = "var1 -> var2", title = "Organelle co-occurrence",
         subtitle = "text: var2 present/absent among cells with var1 vs without var1") +
    theme_minimal()

  p_dumbbell <- ggplot(or_pairs, aes(y = pair)) +
    geom_segment(aes(x = prev_var1_absent, xend = prev_var1_present, yend = pair),
                 colour = "grey70", linewidth = 1) +
    geom_point(aes(x = prev_var1_absent, colour = "without var1"), size = 2.5) +
    geom_point(aes(x = prev_var1_present, colour = "with var1"), size = 2.5) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_colour_manual(values = c("without var1" = "grey50", "with var1" = "#d95f02"), name = NULL) +
    labs(x = "cells with var2", y = "var1 -> var2",
         title = "Frequency of var2 in cells with vs without var1") +
    theme_minimal()

  save_plot(p_forest, "assoc_organelles_odds_ratios", 10, 2 + 0.25 * nrow(or_pairs), opt$plot_dir)
  save_plot(p_dumbbell, "assoc_organelles_frequency_shift", 8, 2 + 0.25 * nrow(or_pairs), opt$plot_dir)
}

## ---- 3) lookup: contingency table of every significant pair ----
save_contingency_tables(to_be_plotted, cor_data, "assoc_organelles_contingency_tables", opt$plot_dir)
