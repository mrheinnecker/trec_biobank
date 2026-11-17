## Plots of the organelle x metadata associations from 07_metadata_associations.R
## 1) heatmap organelles x metadata (all tested pairs, stars = significant)
## 2) composition plot of every significant pair   3) contingency tables of the significant pairs
## Output: plots/assoc_metadata_*

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
results <- read_assoc_results("assoc_metadata", opt$table_dir, "07_metadata_associations.R")
sig <- results %>% filter(p_adj < 0.05)

## ---- 1) heatmap: organelles (rows) x metadata (columns), all tested pairs ----
## Cramer's V for the categorical metadata (0..1), Spearman rho for beauty (-1..1)
n_cells <- map_int(cor_data, ~ if (is_binary(.x)) sum(.x == 1, na.rm = TRUE) else sum(!is.na(.x)))
heat <- results %>%
  filter(!is.na(estimate)) %>%
  mutate(group = feature_group(var1),
         var1 = sprintf("%s (%d)", var1, n_cells[var1]),
         var1 = factor(var1, rev(str_sort(unique(var1), numeric = TRUE))),
         var2 = factor(var2, intersect(meta_feats, unique(var2))))

p_heatmap <- ggplot(heat, aes(var2, var1, fill = estimate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = stars(p_adj)), size = 3, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1),
                       name = "Cramer's V /\nrho (beauty)") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = "organelle (cells with it)", title = "Organelles x metadata",
       subtitle = "* p_adj < 0.05, ** < 0.01, *** < 0.001",
       caption = "Bonferroni within this analysis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        plot.title.position = "plot", strip.text.y = element_text(angle = 0, face = "bold", hjust = 0))

save_plot(p_heatmap, "assoc_metadata_heatmap", 7, 2 + 0.22 * n_distinct(heat$var1), opt$plot_dir)

if (nrow(sig) == 0) {
  message("No significant organelle x metadata pairs")
} else {
  ## ---- 2) composition per metadata category of every significant pair ----
  cat_plots <- sig %>%
    select(var1, var2, p_adj, estimate, test) %>%
    pmap(plot_cat_pair, data = cor_data)
  save_cat_plots(cat_plots, "assoc_metadata_composition", opt$plot_dir)

  ## ---- 3) lookup: contingency tables ----
  save_contingency_tables(sig, cor_data, "assoc_metadata_contingency_tables", opt$plot_dir)
}
