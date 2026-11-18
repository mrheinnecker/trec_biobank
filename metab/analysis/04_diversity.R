# 04: Diversity, separately for 18S eukaryotes and 16S prokaryotes
#   1) alpha diversity (observed ASVs, Shannon) after rarefying every sample to
#      the same read number, per site and size fraction
#   2) beta diversity: Bray-Curtis on relative abundances, PCoA plot
#   3) PERMANOVA: how much of the variation site and size fraction explain
#
# Note: every site was sequenced in exactly one batch (2023 or 2024), so the
# site effect includes any batch effect; the two cannot be separated.

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/trec_biobank/metab/analysis"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages(library(vegan))

asv_info    <- read_prepared("asv_info")
counts_long <- read_prepared("counts_long")
samples     <- read_prepared("samples") %>% filter(pass_qc)

set.seed(1)
min_category_reads <- 2000   # a sample needs this many reads in the category to be included

# Sample x ASV matrix for one category.
category_matrix <- function(cat) {
  m <- counts_long %>%
    filter(run_accession %in% samples$run_accession) %>%
    semi_join(asv_info %>% filter(category == cat), by = "ASV_ID") %>%
    pivot_wider(id_cols = run_accession, names_from = ASV_ID, values_from = reads, values_fill = 0)
  mat <- as.matrix(m[, -1]); rownames(mat) <- m$run_accession
  mat[rowSums(mat) >= min_category_reads, , drop = FALSE]
}

sets <- c(eukaryote_nuclear = "Eukaryotes (18S)", prokaryote = "Prokaryotes (16S)")
alpha_all <- list(); permanova_all <- list()

for (cat in names(sets)) {
  mat <- category_matrix(cat)
  meta <- samples[match(rownames(mat), samples$run_accession), ]
  cat("\n====", sets[[cat]], ":", nrow(mat), "samples,", ncol(mat), "ASVs ====\n")

  # ---- 1. alpha diversity ----------------------------------------------------
  depth <- min(rowSums(mat))
  rare  <- rrarefy(mat, depth)
  alpha <- meta %>%
    select(run_accession, site, date, size_fraction, sequencing_batch) %>%
    mutate(category = cat, rarefied_to = depth,
           observed_asvs = rowSums(rare > 0), shannon = diversity(rare, "shannon"))
  alpha_all[[cat]] <- alpha

  p_alpha <- alpha %>%
    pivot_longer(c(observed_asvs, shannon), names_to = "index") %>%
    ggplot(aes(x = site, y = value, colour = size_fraction)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 0.8) +
    facet_wrap(~ index, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = NULL, colour = "size fraction",
         title = paste(sets[[cat]], "- alpha diversity"), subtitle = paste("rarefied to", depth, "reads")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_alpha, paste0("04_alpha_", cat), width = 10, height = 7)

  # ---- 2. Bray-Curtis + PCoA -------------------------------------------------
  rel <- mat / rowSums(mat)
  bc  <- vegdist(rel, method = "bray")
  pc  <- cmdscale(bc, k = 3, eig = TRUE)
  expl <- round(100 * pc$eig[1:3] / sum(pc$eig[pc$eig > 0]), 1)
  ord <- meta %>% mutate(PC1 = pc$points[, 1], PC2 = pc$points[, 2], PC3 = pc$points[, 3])
  ord %>% select(run_accession, site, date, time_of_day, size_fraction, replicate, PC1, PC2, PC3) %>%
    write_tsv(file.path(table_dir, paste0("04_pcoa_coordinates_", cat, ".tsv")))

  p_pcoa <- ggplot(ord, aes(PC1, PC2, colour = site, shape = size_fraction)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = paste0("PCoA 1 (", expl[1], "%)"), y = paste0("PCoA 2 (", expl[2], "%)"),
         colour = "site", shape = "size fraction",
         title = paste(sets[[cat]], "- Bray-Curtis PCoA"), subtitle = "relative abundances") +
    coord_equal()
  save_plot(p_pcoa, paste0("04_pcoa_", cat), width = 9, height = 7)

  # Same ordination: site as shape, the third axis as colour. The site ->
  # shape mapping is fixed so that it is identical in all plots below.
  site_shapes <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 9, 10, 11, 12, 13, 14)
  all_sites <- sort(unique(samples$site))
  if (length(all_sites) > length(site_shapes)) warning("more sites than shapes; some sites share a shape")
  shape_scale <- scale_shape_manual(values = setNames(rep_len(site_shapes, length(all_sites)), all_sites), drop = FALSE)
  p_pcoa3 <- ggplot(ord, aes(PC1, PC2, shape = site, colour = PC3)) +
    geom_point(size = 2.5, stroke = 0.8) +
    shape_scale +
    scale_colour_viridis_c(option = "plasma", end = 0.95) +
    labs(x = paste0("PCoA 1 (", expl[1], "%)"), y = paste0("PCoA 2 (", expl[2], "%)"),
         colour = paste0("PCoA 3\n(", expl[3], "%)"), shape = "site",
         title = paste(sets[[cat]], "- Bray-Curtis PCoA, axes 1-3"), subtitle = "relative abundances") +
    coord_equal()
  save_plot(p_pcoa3, paste0("04_pcoa3_", cat), width = 9, height = 7)

  # (B) The overall ordination, one panel per size fraction (shared axes).
  p_pcoa3_facet <- ggplot(ord, aes(PC1, PC2, shape = site, colour = PC3)) +
    geom_point(data = select(ord, -size_fraction), colour = "grey88", size = 1.5) +   # all samples in the background
    geom_point(size = 2.5, stroke = 0.8) +
    facet_wrap(~ size_fraction) +
    shape_scale +
    scale_colour_viridis_c(option = "plasma", end = 0.95) +
    labs(x = paste0("PCoA 1 (", expl[1], "%)"), y = paste0("PCoA 2 (", expl[2], "%)"),
         colour = paste0("PCoA 3\n(", expl[3], "%)"), shape = "site",
         title = paste(sets[[cat]], "- Bray-Curtis PCoA of all samples, split by size fraction"),
         subtitle = "same axes in every panel; grey = samples of the other fractions") +
    coord_equal()
  save_plot(p_pcoa3_facet, paste0("04_pcoa3_by_fraction_shared_axes_", cat), width = 14, height = 9)

  # (A) A separate PCoA for each size fraction (own axes, own axis 3 colour).
  min_fraction_samples <- 5
  fraction_plots <- list()
  for (fr in sort(unique(meta$size_fraction))) {
    keep <- meta$size_fraction == fr
    if (sum(keep) < min_fraction_samples) {
      message(sets[[cat]], ": size fraction ", fr, " has only ", sum(keep), " samples, no separate PCoA")
      next
    }
    bc_fr <- vegdist(rel[keep, colSums(rel[keep, , drop = FALSE]) > 0], method = "bray")
    pc_fr <- cmdscale(bc_fr, k = 3, eig = TRUE)
    ex_fr <- round(100 * pc_fr$eig[1:3] / sum(pc_fr$eig[pc_fr$eig > 0]), 1)
    ord_fr <- meta[keep, ] %>% mutate(PC1 = pc_fr$points[, 1], PC2 = pc_fr$points[, 2], PC3 = pc_fr$points[, 3])
    p <- ggplot(ord_fr, aes(PC1, PC2, shape = site, colour = PC3)) +
      geom_point(size = 2.5, stroke = 0.8) +
      shape_scale +
      scale_colour_viridis_c(option = "plasma", end = 0.95) +
      guides(shape = "none", colour = guide_colourbar(barheight = unit(3, "cm"))) +
      labs(x = paste0("PCoA 1 (", ex_fr[1], "%)"), y = paste0("PCoA 2 (", ex_fr[2], "%)"),
           colour = paste0("PCoA 3\n(", ex_fr[3], "%)"),
           title = paste0(fr, " µm  (n = ", sum(keep), ")"))
    # wrap_elements keeps each panel's own colour bar next to it
    fraction_plots[[fr]] <- patchwork::wrap_elements(full = p)
  }
  if (length(fraction_plots) > 0) {
    # One shared site legend for the whole figure (from an otherwise empty plot).
    site_legend <- ggplot(tibble(site = all_sites), aes(0, 0, shape = site)) +
      geom_point(alpha = 0) + shape_scale +
      guides(shape = guide_legend(override.aes = list(alpha = 1, size = 3))) +
      theme_void()
    p_fr <- (patchwork::wrap_plots(fraction_plots, ncol = 2) | site_legend) +
      patchwork::plot_layout(widths = c(1, 0.001), guides = "collect") +
      patchwork::plot_annotation(title = paste(sets[[cat]], "- separate Bray-Curtis PCoA per size fraction"),
                                 subtitle = "each panel has its own axes and its own PCoA 3 colour scale")
    save_plot(p_fr, paste0("04_pcoa3_by_fraction_separate_", cat), width = 14, height = 4.5 * ceiling(length(fraction_plots) / 2))
  }

  # ---- 3. PERMANOVA ----------------------------------------------------------
  # Marginal tests: each term after accounting for the other.
  perm <- adonis2(bc ~ site + size_fraction, data = meta, by = "margin", permutations = 999)
  print(perm)
  permanova_all[[cat]] <- as.data.frame(perm) %>% rownames_to_column("term") %>% mutate(category = cat, .before = 1)
}

write_tsv(bind_rows(alpha_all), file.path(table_dir, "04_alpha_diversity.tsv"))
write_tsv(bind_rows(permanova_all), file.path(table_dir, "04_permanova.tsv"))
cat("plots written to", plot_dir, "\n")
