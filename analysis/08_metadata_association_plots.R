## Plots of the organelle x metadata associations from 07_metadata_associations.R
## 1) heatmap organelles x metadata (all tested pairs, stars = significant)
## 1b) overview in one figure: prevalence per metadata category (direction) + test per organelle x variable
## 2) composition plot of every significant pair   3) contingency tables of the significant pairs
## Output: plots/assoc_metadata_*

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- NA_character_   # NA = read the main sheet from Google Sheets; or a tsv snapshot
default_outdir     <- "C:/projects/trec_biobank"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/trec_biobank/analysis"

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

## ---- 1b) overview: prevalence per metadata category (direction of the associations) ----
## rows = organelles, columns = categories of every metadata variable (one block per variable);
## fill = prevalence in the category minus overall prevalence (percentage points), number = prevalence (%);
## last column of a block ("test") = test of organelle x variable (Cramer's V / rho + stars, Bonferroni)
min_cat_cells <- 5   # categories with fewer cells: prevalence shown in grey (too few cells)
feats_tested <- unique(results$var1)
overview_skip <- c("site")   # sites: own figure (03_site_features.R)
vars_tested  <- setdiff(intersect(meta_feats, unique(results$var2)), overview_skip)
prev <- map_dfr(vars_tested, function(v) {
  d <- cor_data %>% select(cat = all_of(v), all_of(feats_tested)) %>% filter(!is.na(cat)) %>%
    mutate(cat = as.character(cat))
  d %>%
    pivot_longer(-cat, names_to = "feature", values_to = "present") %>%
    filter(!is.na(present)) %>%
    group_by(feature) %>% mutate(overall = mean(present == 1)) %>%
    group_by(feature, cat, overall) %>%
    summarise(n = n(), prev = mean(present == 1), .groups = "drop") %>%
    mutate(variable = v)
}) %>%
  mutate(diff = if_else(n < min_cat_cells, NA_real_, prev - overall))   # too few cells: no colour
cat_n <- prev %>% group_by(variable, cat) %>% summarise(n = max(n), .groups = "drop")
cat_order <- cat_n %>% arrange(variable, str_rank(cat, numeric = TRUE)) %>%
  mutate(lab = sprintf("%s\n(%d)", str_replace(cat, "to", "-"), n))
effect_col <- "test"; p_col <- "p_adj"
fmt_padj <- function(p) case_when(is.na(p) ~ "", p < 0.001 ~ formatC(p, format = "e", digits = 0),
                                  TRUE ~ formatC(p, format = "g", digits = 2))
tests <- results %>% filter(var1 %in% feats_tested, var2 %in% vars_tested) %>%
  transmute(feature = var1, variable = var2, cat = effect_col,
            label = if_else(is.na(estimate), "", paste0(formatC(estimate, digits = 2, format = "f"), stars(p_adj))),
            sig = !is.na(p_adj) & p_adj < 0.05, p_adj)
tests <- bind_rows(tests, tests %>% mutate(cat = p_col, label = fmt_padj(p_adj)))
row_n <- map_int(set_names(feats_tested), ~ sum(cor_data[[.x]] == 1, na.rm = TRUE))
row_levels <- rev(str_sort(sprintf("%s (%d)", feats_tested, row_n[feats_tested]), numeric = TRUE))
shape_rows <- function(d) d %>%
  mutate(group = feature_group(feature),
         row = factor(sprintf("%s (%d)", feature, row_n[feature]), row_levels),
         variable = factor(variable, vars_tested),
         col = factor(paste(variable, cat), c(paste(cat_order$variable, cat_order$cat),
                                              paste(vars_tested, effect_col), paste(vars_tested, p_col))))
prev_plot <- shape_rows(prev)
test_plot <- shape_rows(tests)
col_labels <- c(setNames(rep("test", length(vars_tested)), paste(vars_tested, effect_col)),
                setNames(rep("p_adj", length(vars_tested)), paste(vars_tested, p_col)),
                setNames(cat_order$lab, paste(cat_order$variable, cat_order$cat)))

p_overview <- ggplot(mapping = aes(col, row)) +
  geom_tile(data = prev_plot, aes(fill = diff), colour = "white") +
  geom_text(data = prev_plot, aes(label = round(100 * prev), colour = if_else(n < min_cat_cells, "grey60", "grey10")),
            size = 2.3) +
  geom_tile(data = test_plot, fill = "grey96", colour = "white") +
  geom_text(data = test_plot, aes(label = label, fontface = if_else(sig, "bold", "plain")), size = 2.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-0.6, 0.6),
                       oob = scales::squish, na.value = "grey92", labels = function(x) paste0(ifelse(x > 0, "+", ""), round(100 * x)),
                       name = "prevalence vs.\noverall\n(percentage\npoints)") +
  scale_colour_identity() +
  scale_x_discrete(labels = col_labels) +
  facet_grid(group ~ variable, scales = "free", space = "free", switch = "x") +
  labs(x = NULL, y = "organelle (cells with it)",
       title = "Organelles across metadata categories",
       subtitle = paste0("number = % of the cells of the category with the organelle (grey tile: < ", min_cat_cells,
                         " cells, not coloured); colour = difference to the overall prevalence\n",
                         "column 'test' = Cramer's V (size_frac, sampling_time, life_status) / Spearman rho (beauty); ",
                         "* p_adj < 0.05, ** < 0.01, *** < 0.001;\ncolumn 'p_adj' = p-value (Bonferroni over all organelle x metadata tests of 07, incl. sites); bold = significant; sites: see orgs_per_site"),
       caption = "x labels: category (number of cells)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7, lineheight = 0.85),
        strip.placement = "outside", strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
        plot.title.position = "plot", plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.spacing.x = unit(8, "pt"))

n_cols <- nrow(cat_order) + 2 * length(vars_tested)
save_plot(p_overview, "assoc_metadata_overview", 3 + 0.48 * n_cols, 2.4 + 0.24 * length(feats_tested), opt$plot_dir)
write_tsv(prev %>% select(variable, category = cat, feature, n_cells = n, prevalence = prev, overall, diff),
          file.path(opt$table_dir, "assoc_metadata_prevalence.tsv"))

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
