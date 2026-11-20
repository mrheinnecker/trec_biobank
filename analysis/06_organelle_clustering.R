## Clustering of cells by organelle composition (dinoflagellates only)
## every cell = 0/1 vector of organelles present / absent (same data as the association analysis)
## distance: Jaccard (only shared presences count as similarity, shared absences are ignored:
##           "both lack a pyrenoid" says less than "both have one")
## clustering: hierarchical, Ward linkage (compact clusters of similar size)
##   Ward needs Euclidean distances. Jaccard itself is not Euclidean, but sqrt(Jaccard) is
##   (Gower & Legendre 1986), so Ward runs on sqrt(Jaccard). This is the same as Ward on all
##   PCoA axes of sqrt(Jaccard). With n_axes set, Ward runs on only the first n_axes PCoA axes
##   (drops the small, noisy axes; the distance is otherwise unchanged).
##   Correlated features: they share the leading axes and so weigh more in the distance, with or
##   without n_axes. If that is unwanted, merge / drop redundant features (see clust_feature_correlation).
## number of clusters (k_method, or fixed via k_fixed):
##   "stability": n_boot times a random subset of boot_frac of the cells is clustered again
##     (subsampling without replacement: a bootstrap with replacement would create duplicate cells at
##     distance 0, which hierarchical clustering merges first). For every original cluster the best
##     match in the reclustered subset is found and their overlap measured as Jaccard similarity
##     (Hennig's clusterboot, on the cells in the subset). Mean over the resamples per cluster:
##     >= 0.75 stable, 0.6-0.75 a pattern, < 0.5 dissolves. The agreement of the whole partition is
##     reported as adjusted Rand index (ARI).
##     Chosen k = the largest k in k_range in which every cluster has mean Jaccard >= stab_min
##     (the most detail that still reproduces); if no k qualifies, the k with the highest mean ARI.
##   "silhouette": highest mean silhouette width
## Note: cells of one sample are not independent, so the cluster vs metadata tests are optimistic.

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
## features left out of the clustering only (still in 04/05 and in the cluster profile table):
## pyrenoid and plastoglob are nearly redundant with chloropl (phi 0.72-0.78), so the chloroplast
## counted ~3x in the distance; chloropl stays as the single chloroplast feature
drop_from_clustering <- c("pyrenoid", "plastoglob")
k_range        <- 2:10         # candidate numbers of clusters
k_fixed        <- NA           # set e.g. 4 to override the automatic choice
## env TEM_CLUSTER_K (e.g. from launcher.sh) overrides k_fixed; empty = automatic choice (k_method)
if (Sys.getenv("TEM_CLUSTER_K") != "") k_fixed <- as.integer(Sys.getenv("TEM_CLUSTER_K"))
k_method       <- "stability"  # "stability" or "silhouette" (see header)
n_boot         <- 200          # number of resamples for the stability (env TEM_CLUSTER_BOOT)
if (Sys.getenv("TEM_CLUSTER_BOOT") != "") n_boot <- as.integer(Sys.getenv("TEM_CLUSTER_BOOT"))
boot_frac      <- 0.8          # fraction of cells per resample
stab_min       <- 0.6          # minimum mean Jaccard every cluster needs for a k to count as stable
linkage        <- "ward.D2"     # on sqrt(Jaccard), see header
n_axes         <- NA           # NA = full sqrt(Jaccard) distance; e.g. 5 = Ward on the first 5 PCoA axes
## env TEM_CLUSTER_AXES overrides n_axes
if (Sys.getenv("TEM_CLUSTER_AXES") != "") n_axes <- as.integer(Sys.getenv("TEM_CLUSTER_AXES"))
cor_flag       <- 0.5          # feature pairs with |phi| above this are listed as strongly correlated
meta_cols     <- c("site", "size_frac", "sampling_time", "life_status")

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
  mutate(used = n_present >= min_cells & prevalence <= max_prevalence & !feature %in% drop_from_clustering)
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
d_jac  <- dist(mat, method = "binary")   # = Jaccard distance
d_sqrt <- sqrt(d_jac)                     # Euclidean-embeddable -> valid for Ward
## PCoA of sqrt(Jaccard): all eigenvalues are returned, coordinates for the first max(2, n_axes) axes
pcoa <- cmdscale(d_sqrt, k = max(2, n_axes, na.rm = TRUE), eig = TRUE)
pos_eig <- pcoa$eig[pcoa$eig > 1e-10]
var_expl_all <- pos_eig / sum(pos_eig)
if (!is.na(n_axes)) {
  if (n_axes > ncol(pcoa$points)) stop("n_axes = ", n_axes, " but only ", ncol(pcoa$points), " PCoA axes")
  message(sprintf("clustering on the first %d PCoA axes (%.0f%% of the variation)",
                  n_axes, 100 * sum(var_expl_all[seq_len(n_axes)])))
  d <- dist(pcoa$points[, seq_len(n_axes), drop = FALSE])
} else {
  d <- d_sqrt
}
dist_label <- if (is.na(n_axes)) "sqrt(Jaccard)" else sprintf("first %d PCoA axes of sqrt(Jaccard)", n_axes)
hc <- hclust(d, method = linkage)

sil_by_k <- tibble(k = k_range) %>%
  mutate(mean_sil = map_dbl(k, ~ mean(silhouette(cutree(hc, .x), d)[, "sil_width"])))

## ---- stability of the clusters under resampling (see header) ----
adj_rand <- function(x, y) {
  tab <- table(x, y)
  a <- sum(choose(tab, 2)); b <- sum(choose(rowSums(tab), 2)); c <- sum(choose(colSums(tab), 2))
  e <- b * c / choose(sum(tab), 2)
  (a - e) / ((b + c) / 2 - e)
}
## the distance between two cells does not depend on the other cells -> subset the full matrix
d_mat <- as.matrix(d)
cut_orig <- map(set_names(k_range), ~ cutree(hc, .x))
set.seed(42)
n_sub <- round(boot_frac * nrow(mat))
boot_res <- map_dfr(seq_len(n_boot), function(b) {
  idx <- sort(sample(nrow(mat), n_sub))
  hc_b <- hclust(as.dist(d_mat[idx, idx]), method = linkage)
  map_dfr(k_range, function(kk) {
    orig <- cut_orig[[as.character(kk)]][idx]
    resamp <- cutree(hc_b, kk)
    ## best-matching resampled cluster per original cluster, Jaccard of their cell sets
    jac <- map_dbl(seq_len(kk), function(cl) {
      in_orig <- orig == cl
      if (!any(in_orig)) return(NA_real_)   # cluster not drawn at all in this resample
      max(map_dbl(unique(resamp), ~ sum(in_orig & resamp == .x) / sum(in_orig | resamp == .x)))
    })
    tibble(boot = b, k = kk, cl_raw = seq_len(kk), jaccard = jac, ari = adj_rand(orig, resamp))
  })
})
stab_clusters <- boot_res %>%
  group_by(k, cl_raw) %>%
  summarise(mean_jaccard = mean(jaccard, na.rm = TRUE), n_dissolved = sum(jaccard < 0.5, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(size = map2_int(k, cl_raw, ~ sum(cut_orig[[as.character(.x)]] == .y)))
stab_by_k <- boot_res %>%
  distinct(boot, k, ari) %>%
  group_by(k) %>%
  summarise(mean_ari = mean(ari), ari_low = quantile(ari, 0.025), ari_high = quantile(ari, 0.975)) %>%
  left_join(stab_clusters %>% group_by(k) %>%
              summarise(min_cluster_jaccard = min(mean_jaccard), mean_cluster_jaccard = mean(mean_jaccard)),
            by = "k") %>%
  left_join(sil_by_k, by = "k") %>%
  mutate(stable = min_cluster_jaccard >= stab_min)

k_auto <- if (k_method == "stability") {
  if (any(stab_by_k$stable)) max(stab_by_k$k[stab_by_k$stable]) else stab_by_k$k[which.max(stab_by_k$mean_ari)]
} else {
  sil_by_k$k[which.max(sil_by_k$mean_sil)]
}
k <- if (is.na(k_fixed)) k_auto else k_fixed
k_how <- if (!is.na(k_fixed)) "fixed" else k_method
with(stab_by_k[stab_by_k$k == k, ],
     message(sprintf("k = %d (%s; mean silhouette %.3f, mean ARI %.2f, least stable cluster Jaccard %.2f)",
                     k, k_how, mean_sil, mean_ari, min_cluster_jaccard)))
if (k_method == "stability" && is.na(k_fixed) && !any(stab_by_k$stable))
  message(sprintf("no k with all clusters at mean Jaccard >= %.2f; k chosen by the highest mean ARI", stab_min))

## clusters numbered left to right along the dendrogram, so C1..Ck match the heatmap
cl_raw <- cutree(hc, k)
cl_order <- unique(cl_raw[hc$order])
cells <- cells %>%
  mutate(cluster = factor(paste0("C", match(cl_raw, cl_order)), paste0("C", seq_len(k))),
         sil_width = silhouette(cl_raw, d)[, "sil_width"])

cluster_cols <- setNames(if (k <= 8) scales::brewer_pal(palette = "Dark2")(max(k, 3))[seq_len(k)] else scales::hue_pal()(k),
                         levels(cells$cluster))

## stability tables (cluster names C1..Ck only for the chosen k; other k keep the cutree number)
stab_clusters <- stab_clusters %>%
  mutate(cluster = if_else(k == !!k, paste0("C", match(cl_raw, cl_order)), paste0("k", k, "_", cl_raw)))
write_tsv(stab_by_k, file.path(opt$table_dir, "clust_stability_by_k.tsv"))
write_tsv(stab_clusters %>% select(k, cluster, size, mean_jaccard, n_dissolved) %>% arrange(k, cluster),
          file.path(opt$table_dir, "clust_stability_clusters.tsv"))

## ---- 0) stability plot ----
p_stab_k <- stab_by_k %>%
  select(k, `ARI (whole partition)` = mean_ari, `least stable cluster (Jaccard)` = min_cluster_jaccard,
         `mean over clusters (Jaccard)` = mean_cluster_jaccard) %>%
  pivot_longer(-k, names_to = "measure", values_to = "value") %>%
  ggplot(aes(k, value, colour = measure)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = stab_min, ymax = Inf, fill = "grey92") +
  geom_hline(yintercept = c(0.5, 0.75), linetype = "dotted", colour = "grey50") +
  geom_vline(xintercept = k, colour = "#d95f02", linetype = "dashed") +
  geom_errorbar(data = stab_by_k, aes(k, ymin = ari_low, ymax = ari_high), inherit.aes = FALSE,
                width = 0.15, colour = "#1b9e77") +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("#1b9e77", "#7570b3", "#bdb7e6"), name = NULL) +
  scale_x_continuous(breaks = k_range) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "number of clusters k", y = "agreement with the full clustering",
       title = "Stability of the clustering under resampling",
       subtitle = sprintf("%d resamples of %.0f%% of the cells; ARI with 95%% range\ngrey = every cluster >= %.2f; chosen k = %d (%s)",
                          n_boot, 100 * boot_frac, stab_min, k, k_how)) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "vertical")

p_stab_cl <- boot_res %>%
  filter(k == !!k, !is.na(jaccard)) %>%
  mutate(cluster = factor(paste0("C", match(cl_raw, cl_order)), levels(cells$cluster))) %>%
  ggplot(aes(cluster, jaccard, fill = cluster)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = Inf, fill = "grey92") +
  geom_hline(yintercept = c(0.5, 0.75), linetype = "dotted", colour = "grey50") +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  geom_text(data = stab_clusters %>% filter(k == !!k) %>%
              mutate(cluster = factor(cluster, levels(cells$cluster))),
            aes(cluster, y = 1.04, label = sprintf("%.2f\nn=%d", mean_jaccard, size)),
            inherit.aes = FALSE, size = 2.8, vjust = 0) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Jaccard with best-matching cluster",
       title = sprintf("Per cluster, k = %d", k), subtitle = "mean Jaccard and\ncluster size on top") +
  theme_minimal()

save_plot(p_stab_k + p_stab_cl + plot_layout(widths = c(1.4, 1)), "clust_stability", 12, 5, opt$plot_dir)

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
       title = "Choice of k", subtitle = sprintf("%s, %s linkage; chosen k = %d", dist_label, linkage, k)) +
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

## ---- 3) PCoA of the sqrt(Jaccard) distances (computed in the clustering step) ----
var_expl <- var_expl_all[1:2]
pcoa_df <- cells %>% mutate(PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2])
pcoa_axes <- labs(x = sprintf("PCo1 (%.0f%%)", 100 * var_expl[1]), y = sprintf("PCo2 (%.0f%%)", 100 * var_expl[2]))

p_pcoa_cl <- ggplot(pcoa_df, aes(PCo1, PCo2, colour = cluster)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = cluster_cols) +
  pcoa_axes + labs(title = "PCoA (sqrt Jaccard), coloured by cluster") +
  theme_minimal()
p_pcoa_site <- ggplot(pcoa_df, aes(PCo1, PCo2, colour = site)) +
  geom_jitter(width = 0.01, height = 0.01, size = 1.8, alpha = 0.8) +
  scale_colour_brewer(palette = "Set2") +
  pcoa_axes + labs(title = "coloured by site") +
  theme_minimal()

save_plot(p_pcoa_cl + p_pcoa_site, "clust_pcoa", 12, 5, opt$plot_dir)

## scree: variation per PCoA axis (helps to pick n_axes; the cut is where the curve flattens)
scree <- tibble(axis = seq_along(var_expl_all), var_expl = var_expl_all, cum = cumsum(var_expl_all))
write_tsv(scree, file.path(opt$table_dir, "clust_pcoa_scree.tsv"))
p_scree <- scree %>%
  filter(axis <= 20) %>%
  ggplot(aes(axis, var_expl)) +
  geom_col(aes(fill = !is.na(n_axes) & axis <= n_axes), show.legend = FALSE) +
  geom_line(aes(y = cum / max(cum) * max(var_expl)), colour = "grey40") +
  geom_point(aes(y = cum / max(cum) * max(var_expl)), colour = "grey40", size = 1) +
  scale_fill_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "grey60")) +
  scale_y_continuous(labels = scales::percent,
                     sec.axis = sec_axis(~ . / max(scree$var_expl[scree$axis <= 20]) * max(scree$cum[scree$axis <= 20]),
                                         labels = scales::percent, name = "cumulative")) +
  labs(x = "PCoA axis", y = "variation explained", title = "PCoA of sqrt(Jaccard): scree",
       subtitle = if (is.na(n_axes)) "clustering on the full distance (all axes)"
                  else sprintf("orange = the %d axes used for clustering", n_axes)) +
  theme_minimal()
save_plot(p_scree, "clust_pcoa_scree", 6, 4, opt$plot_dir)

## ---- correlated features (phi = Pearson on 0/1) ----
## strongly correlated features carry the same information twice and weigh more in the clustering
phi <- cor(mat)
phi_pairs <- as_tibble(as.table(phi), .name_repair = ~ c("feature1", "feature2", "phi")) %>%
  filter(as.character(feature1) < as.character(feature2)) %>%
  arrange(desc(abs(phi)))
write_tsv(phi_pairs, file.path(opt$table_dir, "clust_feature_correlation.tsv"))
strong <- phi_pairs %>% filter(abs(phi) > cor_flag)
if (nrow(strong) > 0)
  message(sprintf("%d feature pairs with |phi| > %.2f: %s", nrow(strong), cor_flag,
                  paste0(strong$feature1, "~", strong$feature2, " (", round(strong$phi, 2), ")", collapse = ", ")))

feat_ord <- colnames(mat)[hclust(as.dist(1 - abs(phi)), method = "average")$order]
p_phi <- as_tibble(as.table(phi), .name_repair = ~ c("feature1", "feature2", "phi")) %>%
  mutate(across(c(feature1, feature2), ~ factor(.x, feat_ord))) %>%
  ggplot(aes(feature1, feature2, fill = phi)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = if_else(abs(phi) > cor_flag & feature1 != feature2, sprintf("%.2f", phi), "")), size = 2.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "phi") +
  labs(x = NULL, y = NULL, title = "Correlation between the clustering features",
       subtitle = sprintf("phi coefficient; values shown for |phi| > %.2f", cor_flag)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_plot(p_phi, "clust_feature_correlation", 3 + 0.3 * ncol(mat), 2.5 + 0.3 * ncol(mat), opt$plot_dir)

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

save_plot(p_profile, "clust_profiles", 3.5 + 0.8 * k, 2 + 0.22 * length(used_feats), opt$plot_dir)

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
