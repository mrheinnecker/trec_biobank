## Clustering of cells by organelle composition (dinoflagellates only)
## every cell = 0/1 vector of organelles present / absent (same data as the association analysis)
## distance: Jaccard (only shared presences count as similarity, shared absences are ignored:
##           "both lack a pyrenoid" says less than "both have one")
## clustering: hierarchical, Ward linkage (compact clusters of similar size)
## number of clusters: highest mean silhouette width over k_range, or fixed via k_fixed
## Note: cells of one sample are not independent, so the cluster vs metadata tests are optimistic.

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- "C:/projects/tem_classification/all_annotations_2026-10-15.tsv"
default_outdir     <- "C:/projects/tem_classification"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/tem_classification/analysis"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "assoc_common.R"))
suppressPackageStartupMessages({
  library(cluster)
  library(ComplexHeatmap)
  library(patchwork)
})
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)

## ---- settings ----
## organelles to cluster on: remove / comment out what should not be used
## (features that are too rare / too common are dropped afterwards via min_cells / max_prevalence)
use_features <- c(
  #"cell_cover", "nucleus", "nucleolus", 
  "tricho", "core_vesic", 
  #"flagell_app", "symbiosis",
                  "chloropl", "plastoglob", 
  #"pyrenoid", "pusule", 
  "retic_net", "electr_sheets", "tubul_net",
                  "food_vac", "starch", "cryst_rich_RN", 
  #"eyespot", 
  "rhabdo", #"ER", "mito", "golgi",
                  "lipid_drop", 
  #"fibrous_body", "put_division", "SDV", "lysosome", 
  "SER_whorls"
  )
include_undescribed <- TRUE    # also cluster on the undescribed_<n> structures
min_cells      <- 10            # drop features present in fewer cells (no information, only noise)
max_prevalence <- 0.95         # drop features present in more than this fraction of cells (mito, ...)
k_range        <- 2:10         # candidate numbers of clusters
k_fixed        <- NA           # set e.g. 4 to override the silhouette choice
linkage        <- "ward.D2"
meta_cols      <- c("site", "size_frac", "sampling_time", "life_status")

## ---- data ----
dino <- load_dino_data(opt$input)
cor_data <- build_cor_data(dino)

## organelle columns (retic_net already corrected for cryst_rich_RN in load_dino_data)
unknown <- setdiff(use_features, names(dino))
if (length(unknown) > 0) stop("use_features not in the data: ", paste(unknown, collapse = ", "))
feature_cols <- use_features
if (include_undescribed) feature_cols <- c(feature_cols, undescribed_cols(dino))

cells <- cor_data %>%
  mutate(name = dino$name) %>%
  select(name, all_of(meta_cols), all_of(feature_cols))

n_before <- nrow(cells)
cells <- cells %>% drop_na(all_of(feature_cols))
message(sprintf("%d cells, %d dropped because of missing feature values", nrow(cells), n_before - nrow(cells)))

prevalence <- colMeans(cells[feature_cols])
feature_info <- tibble(feature = feature_cols, n_present = colSums(cells[feature_cols]), prevalence = prevalence) %>%
  mutate(used = n_present >= min_cells & prevalence <= max_prevalence)
used_feats <- feature_info %>% filter(used) %>% pull(feature)
message("features not used: ", paste(setdiff(feature_cols, used_feats), collapse = ", "))

mat <- cells %>% select(all_of(used_feats)) %>% as.matrix()
rownames(mat) <- cells$name

## cells without any of the used features cannot be placed by Jaccard (distance undefined)
empty <- rowSums(mat) == 0
if (any(empty)) {
  message(sprintf("%d cells without any used feature dropped", sum(empty)))
  mat <- mat[!empty, , drop = FALSE]
  cells <- cells[!empty, ]
}

## ---- clustering ----
d  <- dist(mat, method = "binary")   # = Jaccard distance
hc <- hclust(d, method = linkage)

sil_by_k <- tibble(k = k_range) %>%
  mutate(mean_sil = map_dbl(k, ~ mean(silhouette(cutree(hc, .x), d)[, "sil_width"])))
k <- if (is.na(k_fixed)) sil_by_k$k[which.max(sil_by_k$mean_sil)] else k_fixed
message(sprintf("k = %d (mean silhouette %.3f)", k, sil_by_k$mean_sil[sil_by_k$k == k]))

## clusters numbered left to right along the dendrogram, so C1..Ck match the heatmap
cl_raw <- cutree(hc, k)
cl_order <- unique(cl_raw[hc$order])
cells <- cells %>%
  mutate(cluster = factor(paste0("C", match(cl_raw, cl_order)), paste0("C", seq_len(k))),
         sil_width = silhouette(cl_raw, d)[, "sil_width"])

cluster_cols <- setNames(if (k <= 8) scales::brewer_pal(palette = "Dark2")(max(k, 3))[seq_len(k)] else scales::hue_pal()(k),
                         levels(cells$cluster))

## ---- tables ----
## fraction of cells with each feature per cluster (all features, also the unused ones)
cluster_profiles <- cells %>%
  pivot_longer(all_of(feature_cols), names_to = "feature", values_to = "present") %>%
  group_by(cluster, feature) %>%
  summarise(n_cells = n(), n_present = sum(present), frac = mean(present), .groups = "drop") %>%
  left_join(feature_info %>% select(feature, used, overall = prevalence), by = "feature")

## cluster vs metadata: fisher test with simulated p (few tests -> BH)
set.seed(42)
meta_tests <- map_dfr(meta_cols, function(v) {
  tab <- table(cells[[v]], cells$cluster)
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  tibble(variable = v, n = sum(tab),
         cramers_v = unname(sqrt(chi2 / (sum(tab) * (min(dim(tab)) - 1)))),
         p_value = fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value)
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

write_tsv(cells %>% select(name, all_of(meta_cols), cluster, sil_width), file.path(opt$table_dir, "clust_assignments.tsv"))
write_tsv(cluster_profiles, file.path(opt$table_dir, "clust_profiles.tsv"))
write_tsv(sil_by_k, file.path(opt$table_dir, "clust_silhouette_by_k.tsv"))
write_tsv(feature_info, file.path(opt$table_dir, "clust_features.tsv"))
write_tsv(meta_tests, file.path(opt$table_dir, "clust_metadata_tests.tsv"))

## ---- 1) choice of k ----
p_sil <- ggplot(sil_by_k, aes(k, mean_sil)) +
  geom_line(colour = "grey50") +
  geom_point(aes(colour = k == !!k), size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "grey30")) +
  scale_x_continuous(breaks = k_range) +
  labs(x = "number of clusters k", y = "mean silhouette width",
       title = "Choice of k", subtitle = sprintf("Jaccard distance, %s linkage; chosen k = %d", linkage, k)) +
  theme_minimal()

## silhouette of every cell in the chosen clustering: negative = closer to another cluster
p_sil_cells <- cells %>%
  arrange(cluster, desc(sil_width)) %>%
  mutate(idx = row_number()) %>%
  ggplot(aes(idx, sil_width, fill = cluster)) +
  geom_col(width = 1) +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "cells (sorted by cluster and silhouette)", y = "silhouette width",
       title = "Silhouette per cell") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())

save_plot(p_sil + p_sil_cells + plot_layout(widths = c(1, 1.6)), "clust_silhouette", 11, 4, opt$plot_dir)

## ---- 2) heatmap: cells x features, dendrogram, cluster split, metadata on top ----
## features in the same groups as the other analyses: chloroplast, common, specific
feat_group <- case_when(used_feats %in% chloro_feats ~ "chloroplast",
                        used_feats %in% common_feats ~ "common",
                        str_detect(used_feats, "^undescribed_") ~ "undescribed",
                        TRUE ~ "specific") %>%
  factor(c("chloroplast", "common", "specific", "undescribed"))

meta_palette <- function(x, pal) {
  lv <- sort(unique(na.omit(x)))
  setNames(scales::brewer_pal(palette = pal)(max(length(lv), 3))[seq_along(lv)], lv)
}
top_anno <- HeatmapAnnotation(
  cluster       = cells$cluster,
  site          = cells$site,
  size_frac     = cells$size_frac,
  sampling_time = cells$sampling_time,
  life_status   = cells$life_status,
  col = list(cluster       = cluster_cols,
             site          = meta_palette(cells$site, "Set2"),
             size_frac     = meta_palette(cells$size_frac, "Pastel1"),
             sampling_time = meta_palette(cells$sampling_time, "Accent"),
             life_status   = meta_palette(cells$life_status, "Greys")),
  na_col = "white", annotation_name_side = "left", simple_anno_size = unit(3, "mm")
)
## prevalence of each feature on the right
right_anno <- rowAnnotation(
  `% cells` = anno_barplot(colMeans(mat) * 100, width = unit(2, "cm"), gp = gpar(fill = "grey60", col = NA))
)

ht <- Heatmap(t(mat), name = "organelle",
              col = c(`0` = "grey95", `1` = "#2b2b2b"),
              heatmap_legend_param = list(at = c(0, 1), labels = c("absent", "present")),
              cluster_columns = as.dendrogram(hc), column_split = k,
              column_title = levels(cells$cluster), column_gap = unit(1, "mm"),
              column_dend_height = unit(2, "cm"),
              cluster_rows = TRUE, clustering_distance_rows = "binary", row_split = feat_group,
              row_title_rot = 0, row_names_gp = gpar(fontsize = 8),
              show_column_names = FALSE,
              top_annotation = top_anno, right_annotation = right_anno)

for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf(file.path(opt$plot_dir, "clust_heatmap.pdf"), width = 14, height = 7)
  else png(file.path(opt$plot_dir, "clust_heatmap.png"), width = 14, height = 7, units = "in", res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
}

## ---- 3) PCoA of the Jaccard distances ----
pcoa <- cmdscale(d, k = 2, eig = TRUE)
var_expl <- pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0])
pcoa_df <- cells %>% mutate(PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2])
pcoa_axes <- labs(x = sprintf("PCo1 (%.0f%%)", 100 * var_expl[1]), y = sprintf("PCo2 (%.0f%%)", 100 * var_expl[2]))

p_pcoa_cl <- ggplot(pcoa_df, aes(PCo1, PCo2, colour = cluster)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = cluster_cols) +
  pcoa_axes + labs(title = "PCoA (Jaccard), coloured by cluster") +
  theme_minimal()
p_pcoa_site <- ggplot(pcoa_df, aes(PCo1, PCo2, colour = site)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1.8, alpha = 0.8) +
  scale_colour_brewer(palette = "Set2") +
  pcoa_axes + labs(title = "coloured by site") +
  theme_minimal()

save_plot(p_pcoa_cl + p_pcoa_site, "clust_pcoa", 12, 5, opt$plot_dir)

## ---- 4) cluster profiles: what defines each cluster ----
profile_plot_data <- cluster_profiles %>%
  filter(used) %>%
  mutate(feature = factor(feature, rev(used_feats[order(feat_group)])),
         group = feat_group[match(as.character(feature), used_feats)])
n_per_cluster <- count(cells, cluster)

p_profile <- ggplot(profile_plot_data, aes(cluster, feature, fill = frac)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f", 100 * frac), colour = frac > 0.6), size = 2.8) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_viridis_c(option = "mako", direction = -1, labels = scales::percent, limits = c(0, 1),
                       name = "cells with\nfeature") +
  scale_x_discrete(labels = setNames(sprintf("%s\n(n=%d)", n_per_cluster$cluster, n_per_cluster$n), n_per_cluster$cluster)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "Organelle profile per cluster", subtitle = "% of cells in the cluster with the feature") +
  theme_minimal() +
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0, face = "bold"))

save_plot(p_profile, "clust_profiles", 2.5 + 0.7 * k, 2 + 0.22 * length(used_feats), opt$plot_dir)

## ---- 5) clusters vs metadata ----
meta_long <- cells %>%
  select(cluster, all_of(meta_cols)) %>%
  pivot_longer(all_of(meta_cols), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  count(variable, value, cluster) %>%
  left_join(meta_tests %>% select(variable, cramers_v, p_adj), by = "variable") %>%
  mutate(variable = factor(sprintf("%s\nV = %.2f, p_adj = %.2g", variable, cramers_v, p_adj),
                           unique(sprintf("%s\nV = %.2f, p_adj = %.2g", variable, cramers_v, p_adj))))

p_meta <- ggplot(meta_long, aes(value, n, fill = cluster)) +
  geom_col(position = "fill") +
  geom_text(data = meta_long %>% group_by(variable, value) %>% summarise(n = sum(n), .groups = "drop"),
            aes(value, y = 1.02, label = n), inherit.aes = FALSE, size = 2.8, vjust = 0) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "fraction of cells", title = "Cluster composition per metadata category",
       subtitle = "Cramer's V + fisher test (simulated p, BH); n cells on top; cells treated as independent") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_meta, "clust_metadata", 12, 4.5, opt$plot_dir)
