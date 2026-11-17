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
default_input      <- NA_character_   # NA = read the main sheet from Google Sheets; or a tsv snapshot
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
## the organelles to cluster on (use_features, include_undescribed) are set in assoc_common.R,
## shared with the association analysis (04/05)
## (features that are too rare / too common are dropped afterwards via min_cells / max_prevalence)
min_cells      <- 10            # drop features present in fewer cells (no information, only noise)
max_prevalence <- 0.95         # drop features present in more than this fraction of cells (mito, ...)
k_range        <- 2:10         # candidate numbers of clusters
k_fixed        <- NA           # set e.g. 4 to override the silhouette choice
## env TEM_CLUSTER_K (e.g. from launcher.sh) overrides k_fixed; empty = silhouette choice
if (Sys.getenv("TEM_CLUSTER_K") != "") k_fixed <- as.integer(Sys.getenv("TEM_CLUSTER_K"))
linkage        <- "ward.D2"
meta_cols      <- c("site", "size_frac", "sampling_time", "life_status")

## ---- data ----
dino <- load_dino_data(opt)
cor_data <- build_cor_data(dino)

## organelle columns (retic_net already corrected for cryst_rich_RN in load_dino_data)
## (features in < min_feature_cells cells are already removed by assoc_features)
feature_cols <- assoc_features(dino, verbose = TRUE)

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

## MoBIE collection table: one view/grid per cluster (overview figure: step 10);
## also written as tab clusters_k<k> into the collection spreadsheet (TEM_COLLECTION_SHEET_URL)
cluster_collection <- dino %>%
  select(any_of(collection_cols)) %>%
  inner_join(cells %>% select(name, cluster, sil_width), by = "name") %>%
  arrange(cluster, site, name) %>%
  mutate(group = as.character(cluster), view = group, grid = group, sil_width = round(sil_width, 3)) %>%
  select(any_of(collection_cols), group, sil_width) %>%
  left_join(cells %>% select(name, all_of(used_feats)), by = "name")
collection_dir <- file.path(opt$table_dir, "collections")
dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
cluster_tab <- paste0("clusters_k", k)
write_tsv(cluster_collection, file.path(collection_dir, paste0(cluster_tab, ".tsv")), na = "")
push_collection_tabs(set_names(list(cluster_collection), cluster_tab), opt)

## tree of the clustering for the combined tree + images figure of step 10:
## leaves (cell, position in the tree order, cluster, colour) and merges (hclust convention:
## negative = leaf number, positive = earlier merge step)
write_tsv(tibble(name = hc$labels, leaf_pos = match(seq_along(hc$labels), hc$order),
                 cluster = as.character(cells$cluster), color = unname(cluster_cols[as.character(cells$cluster)])),
          file.path(opt$table_dir, paste0("clust_tree_k", k, "_leaves.tsv")))
write_tsv(tibble(step = seq_len(nrow(hc$merge)), left = hc$merge[, 1], right = hc$merge[, 2], height = hc$height),
          file.path(opt$table_dir, paste0("clust_tree_k", k, "_merges.tsv")))

## ---- taxonomic class (taxo_class, filled for some cells) per cluster ----
## cells of the analysis with a taxo_class; cells dropped before clustering -> "not clustered"
taxo_cells <- dino %>%
  filter(!is.na(taxo_class), str_trim(as.character(taxo_class)) != "") %>%
  select(name, site, taxo_class) %>%
  left_join(cells %>% select(name, cluster), by = "name") %>%
  mutate(cluster = fct_na_value_to_level(factor(cluster, levels(cells$cluster)), "not clustered")) %>%
  arrange(taxo_class, cluster, name)
write_tsv(taxo_cells, file.path(opt$table_dir, "clust_taxo_class_cells.tsv"))

if (nrow(taxo_cells) > 0) {
  taxo_counts <- taxo_cells %>%
    count(taxo_class, cluster, .drop = FALSE) %>%
    group_by(taxo_class) %>%
    mutate(frac = n / sum(n), taxo_label = sprintf("%s (n = %d)", taxo_class, sum(n))) %>%
    ungroup()
  write_tsv(taxo_counts %>% select(taxo_class, cluster, n) %>% pivot_wider(names_from = cluster, values_from = n),
            file.path(opt$table_dir, "clust_taxo_class_counts.tsv"))

  p_taxo <- ggplot(taxo_counts, aes(cluster, fct_rev(taxo_label), fill = frac)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = if_else(n > 0, as.character(n), ""), colour = frac > 0.6), size = 3.5) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1a1a1a"), guide = "none") +
    scale_fill_gradient(low = "#f4f8fd", high = "#184f95", limits = c(0, 1), labels = scales::percent,
                        name = "share of\nthe class") +
    labs(x = NULL, y = NULL, title = "Taxonomic class per cluster",
         subtitle = sprintf("dinoflagellates with taxo_class, k = %d | numbers = cells", k)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), plot.title.position = "plot")

  save_plot(p_taxo, "clust_taxo_class", 2.5 + 0.7 * nlevels(taxo_counts$cluster),
            1.5 + 0.4 * n_distinct(taxo_counts$taxo_class), opt$plot_dir)
}

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

## the same annotation tracks + presence/absence for the combined tree figure of step 10:
## per cell the metadata values with their colours (as in the heatmap) and the 0/1 organelles;
## per organelle its group, row order (clustered within the group, as in the heatmap) and prevalence
track_palettes <- list(cluster       = cluster_cols,
                       site          = meta_palette(cells$site, "Set2"),
                       size_frac     = meta_palette(cells$size_frac, "Pastel1"),
                       sampling_time = meta_palette(cells$sampling_time, "Accent"),
                       life_status   = meta_palette(cells$life_status, "Greys"))
tracks <- tibble(name = cells$name)
for (tr in names(track_palettes)) {
  v <- as.character(cells[[tr]])
  tracks[[tr]] <- v
  tracks[[paste0(tr, "_color")]] <- unname(coalesce(track_palettes[[tr]][v], "#FFFFFF"))
}
## taxonomic class where annotated (bottom track + coloured frame around the image in step 10)
taxo <- dino %>% select(name, taxo_class) %>%
  mutate(taxo_class = na_if(str_trim(as.character(taxo_class)), ""))
taxo_lv <- sort(unique(na.omit(taxo$taxo_class)))
taxo_pal <- setNames(c("#E41A1C", "#377EB8", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#4DAF4A")[seq_along(taxo_lv)], taxo_lv)
tracks <- tracks %>%
  left_join(taxo, by = "name") %>%
  mutate(taxo_class_color = unname(coalesce(taxo_pal[taxo_class], "#FFFFFF")))
tracks <- bind_cols(tracks, as_tibble(mat))
write_tsv(tracks, file.path(opt$table_dir, paste0("clust_tree_k", k, "_tracks.tsv")), na = "")

row_info <- tibble(feature = used_feats, group = feat_group, prevalence = colMeans(mat)) %>%
  group_by(group) %>%
  mutate(order_in_group = if (n() > 2) {
    match(feature, feature[hclust(dist(t(mat[, feature, drop = FALSE]), method = "binary"))$order])
  } else seq_len(n())) %>%
  ungroup() %>%
  arrange(group, order_in_group)
write_tsv(row_info, file.path(opt$table_dir, paste0("clust_tree_k", k, "_rows.tsv")))

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
