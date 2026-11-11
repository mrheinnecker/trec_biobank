## Features across sampling sites, with site tests (dinoflagellates only)
## site_feature_plot(): heatmap of the cells with each 0/1 feature per site + p-values of three site tests
## NAs dropped per feature, BH correction over the features (per test)
## KW:   Kruskal-Wallis on the 0/1 values; for 0/1 data it is the chi-squared test of the
##       site x present table (H = (N-1)/N * X^2), so KW and chi2 give almost the same p
## chi2: chi-squared test of the site x present table
##       both treat every cell as independent -> optimistic p-values
## GLMM: binomial mixed model present ~ site + (1 | sample), likelihood-ratio test against
##       present ~ 1 + (1 | sample); the sample (barcode) takes up the cells of one sample being alike,
##       so site is tested against the sample-to-sample variation
## outlines: for features with a significant GLMM site test, one model per site,
##       present ~ this site vs. all others + (1 | sample), Wald test of the site term
##       p < 0.05 -> that site has clearly more (solid) / fewer (dashed) cells with the feature

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(patchwork)
})
opt <- parse_args()

dino <- load_dino_data(opt$input)

glmm_site_p <- function(d) {
  ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  tryCatch({
    m1 <- suppressMessages(suppressWarnings(glmer(value ~ site + (1 | sample), data = d, family = binomial, control = ctrl)))
    m0 <- suppressMessages(suppressWarnings(glmer(value ~ 1 + (1 | sample), data = d, family = binomial, control = ctrl)))
    anova(m0, m1)$`Pr(>Chisq)`[2]
  }, error = function(e) NA_real_)
}

## feat_groups: named list group -> features, top to bottom; features not listed go into "specific"
## (numeric order). One group only (e.g. undescribed) -> no group labels
site_feature_plot <- function(annot, feature_cols, feat_groups = list()) {
  tot_cells_per_site <- annot %>%
    group_by(site) %>%
    tally()

  feat_long <- annot %>%
    mutate(sample = str_extract(source_name, "(?<=_)\\d{9}(?=_[A-Z]{3}_)")) %>%
    select(site, sample, all_of(feature_cols)) %>%
    pivot_longer(cols = all_of(feature_cols)) %>%
    filter(!is.na(value))
  ## row order + group of every feature, the same in all panels
  ## numbers in the names sorted numerically (undescribed_2 before undescribed_11)
  feat_groups <- map(feat_groups, ~ intersect(.x, feature_cols))
  feat_groups$specific <- str_sort(setdiff(feature_cols, unlist(feat_groups)), numeric = TRUE)
  feat_groups <- compact(feat_groups)
  feat_info <- enframe(feat_groups, name = "group", value = "name") %>% unnest(name)
  add_feat_info <- function(d) {
    d %>%
      mutate(name = as.character(name)) %>%
      left_join(feat_info, by = "name") %>%
      mutate(name  = factor(name, rev(feat_info$name)),  # first feature on top
             group = factor(group, names(feat_groups)))
  }
  ## one block per group, the group name only on the right panel
  group_facets <- facet_grid(group ~ ., scales = "free_y", space = "free_y")
  no_strip <- theme(strip.text.y = element_blank(), strip.background = element_blank())
  if (length(feat_groups) == 1) group_facets <- NULL

  across_sites <- feat_long %>%
    group_by(site, name) %>%
    summarize(n_abs = sum(value), .groups = "drop") %>%
    left_join(tot_cells_per_site, by = "site") %>%
    mutate(frac = n_abs / n)

  site_tests <- feat_long %>%
    group_by(name) %>%
    filter(n_distinct(value) > 1) %>%
    summarise(
      KW   = kruskal.test(value ~ factor(site))$p.value,
      chi2 = suppressWarnings(chisq.test(table(value, site))$p.value),
      GLMM = glmm_site_p(cur_data()),
      .groups = "drop"
    ) %>%
    pivot_longer(c(KW, chi2, GLMM), names_to = "test", values_to = "p") %>%
    group_by(test) %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    mutate(test = factor(test, levels = c("KW", "chi2", "GLMM")),
           p_class = cut(p_adj, c(0, 0.001, 0.01, 0.05, 1), include.lowest = TRUE,
                         labels = c("< 0.001", "< 0.01", "< 0.05", "n.s."))) %>%
    ## features without a test (never present) keep an empty row
    complete(name = feature_cols, test)

  glmm_sign <- site_tests %>% filter(test == "GLMM", p_adj < 0.05) %>% pull(name)

  glmm_site_effects <- feat_long %>%
    filter(name %in% glmm_sign) %>%
    group_by(name) %>%
    group_modify(~ {
      d <- .x
      map_dfr(sort(unique(d$site)), function(s) {
        d$this_site <- as.integer(d$site == s)
        tryCatch({
          m <- suppressMessages(suppressWarnings(glmer(value ~ this_site + (1 | sample), data = d, family = binomial,
                                                       control = glmerControl(optimizer = "bobyqa"))))
          co <- summary(m)$coefficients["this_site", ]
          tibble(site = s, glmm_est = co[["Estimate"]], glmm_p = co[["Pr(>|z|)"]])
        }, error = function(e) tibble(site = s, glmm_est = NA_real_, glmm_p = NA_real_))
      })
    }) %>%
    ungroup()
  if (nrow(glmm_site_effects) == 0) {
    glmm_site_effects <- tibble(name = character(), site = character(), glmm_est = numeric(), glmm_p = numeric())
  }
  glmm_site_effects <- glmm_site_effects %>%
    mutate(glmm_mark = case_when(glmm_p < 0.05 & glmm_est > 0 ~ "more",
                                 glmm_p < 0.05 & glmm_est < 0 ~ "fewer"))

  across_sites <- across_sites %>%
    left_join(glmm_site_effects %>% select(name, site, glmm_mark), by = c("name", "site"))

  across_sites <- add_feat_info(across_sites)

  p_feat <- ggplot(across_sites, aes(x = site, y = name, fill = frac)) +
    geom_tile() +
    geom_tile(data = filter(across_sites, !is.na(glmm_mark)),
              aes(linetype = glmm_mark), fill = NA, colour = "#00E5FF", linewidth = 0.8,
              width = 0.92, height = 0.88) +
    ## black text on the bright end of the colour scale (relative to the max of this plot)
    geom_text(aes(label = n_abs, colour = ifelse(frac > 0.75 * max(frac), "black", "white"))) +
    scale_fill_viridis_c(option = "inferno") +
    scale_colour_identity() +
    scale_linetype_manual(values = c(more = "solid", fewer = "22"), drop = FALSE,
                          limits = c("more", "fewer"), name = "site effect (GLMM sign.)") +
    group_facets + no_strip

  ## p-values of the three site tests, BH-adjusted, colour = significance class
  p_pval <- ggplot(add_feat_info(site_tests), aes(x = test, y = name, fill = p_class)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = case_when(is.na(p_adj) ~ "",
                                    p_adj < 0.001 ~ sprintf("%.0e", p_adj),
                                    TRUE ~ formatC(p_adj, digits = 1, format = "g")),
                  colour = ifelse(p_class %in% c("< 0.001", "< 0.01"), "white", "black")), size = 2.1) +
    scale_colour_identity() +
    scale_fill_manual(values = c(`< 0.001` = "#08519C", `< 0.01` = "#4292C6", `< 0.05` = "#9ECAE1", n.s. = "grey92"),
                      drop = FALSE, na.translate = FALSE, name = "site test\n(BH-adj. p)") +
    labs(x = NULL) +
    group_facets + no_strip +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )

  tot_per_feature <- across_sites %>%
    group_by(name) %>%
    summarize(n_tot = sum(n_abs)) %>%
    add_feat_info()

  p_bar <- ggplot(tot_cells_per_site, aes(y = n, x = site)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.3, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )

  p_side <- ggplot(tot_per_feature, aes(x = n_tot, y = name)) +
    geom_col() +
    geom_text(aes(label = n_tot), hjust = -0.2, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    group_facets +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )

  p_comb <- p_bar + plot_spacer() + plot_spacer() + p_feat + p_pval + p_side +
    plot_layout(ncol = 3, widths = c(3, 1.1, 1), heights = c(0.5, 2), guides = "collect")

  list(plot = p_comb, across_sites = across_sites, site_tests = site_tests, glmm_site_effects = glmm_site_effects)
}

## writes the plot + the three result tables under one prefix
save_site_feature_results <- function(res, prefix, width, height) {
  save_plot(res$plot, prefix, width, height, opt$plot_dir)
  write_tsv(res$across_sites, file.path(opt$table_dir, paste0(prefix, "_counts.tsv")))
  write_tsv(res$site_tests, file.path(opt$table_dir, paste0(prefix, "_site_tests.tsv")))
  write_tsv(res$glmm_site_effects, file.path(opt$table_dir, paste0(prefix, "_glmm_site_effects.tsv")))
}

## organelles: chloroplast related, common organelles, the rest specific
## (same groups as the association heatmap)
org_res <- site_feature_plot(dino, organelle_cols(dino),
                             feat_groups = list(chloroplast = chloro_feats, common = common_feats))
save_site_feature_results(org_res, "orgs_per_site", 9, 9)

## undescribed organelles
udo_cols <- undescribed_cols(dino)
udo_res <- site_feature_plot(dino, udo_cols)
save_site_feature_results(udo_res, "udos_per_site", 9, 3 + 0.25 * length(udo_cols))
