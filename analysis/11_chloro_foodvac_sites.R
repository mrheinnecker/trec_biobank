## Trophic groups across sites (dinoflagellates only)
## every cell with chloropl and food_vac annotated goes into one of three groups:
##   chloroplast + food vacuole   (chloropl = 1, food_vac = 1)
##   chloroplast only             (chloropl = 1, food_vac = 0)
##   no chloroplast               (chloropl = 0, with or without food vacuole)
## cells with a missing chloropl / food_vac value are left out (uncertain "?" already counts as 0)
## Plots:  A) fraction of each group per site (n cells in the segments, samples / dates under the site)
##         B) the same per sample (barcode), to see if a site pattern holds in all its samples
## Tests:  site x group: fisher test with simulated p (every cell independent -> optimistic), Cramer's V
##         per group: binomial GLMM  in_group ~ site + (1 | sample)  vs  in_group ~ 1 + (1 | sample),
##         likelihood-ratio test, BH over the three groups (takes up that cells of one sample are alike)
## Output: plots/chloro_foodvac_sites, tables/chloro_foodvac_{cells,site_counts,sample_counts,tests}.tsv

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- NA_character_   # NA = read the main sheet from Google Sheets; or a tsv snapshot
default_outdir     <- "C:/projects/trec_biobank"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/trec_biobank/analysis"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(patchwork)
})
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)

## ---- settings ----
group_levels <- c("chloroplast + food vacuole", "chloroplast only", "no chloroplast")
group_cols   <- c(`chloroplast + food vacuole` = "#D95F02",
                  `chloroplast only`           = "#1B9E77",
                  `no chloroplast`             = "#7570B3")

## ---- data ----
dino <- load_dino_data(opt) %>%
  mutate(sample        = str_extract(source_name, "(?<=_)\\d{9}(?=_[A-Z]{3}_)"),
         sampling_date = str_match(source_name, "_\\d{9}_[A-Z]{3}_[^_]+_(\\d{8})_")[, 2],
         sampling_date = format(as.Date(sampling_date, "%Y%m%d"), "%Y-%m-%d"))

n_all <- nrow(dino)
cells <- dino %>%
  filter(!is.na(chloropl), !is.na(food_vac)) %>%
  mutate(troph = case_when(chloropl == 1 & food_vac == 1 ~ group_levels[1],
                           chloropl == 1                 ~ group_levels[2],
                           TRUE                          ~ group_levels[3]),
         troph = factor(troph, group_levels)) %>%
  select(name, site, sample, sampling_date, size_frac, sampling_time, chloropl, food_vac, troph)
message(sprintf("%d cells, %d without chloropl / food_vac value left out", nrow(cells), n_all - nrow(cells)))

## ---- counts ----
site_counts <- cells %>%
  count(site, troph, .drop = FALSE) %>%
  group_by(site) %>%
  mutate(n_site = sum(n), frac = n / n_site) %>%
  ungroup()

sample_counts <- cells %>%
  count(site, sample, sampling_date, size_frac, sampling_time, troph, .drop = FALSE) %>%
  group_by(sample) %>%
  mutate(n_sample = sum(n), frac = n / n_sample) %>%
  ungroup() %>%
  filter(n_sample > 0)

samples_per_site <- cells %>%
  group_by(site) %>%
  summarise(n_samples = n_distinct(sample), n_dates = n_distinct(sampling_date), .groups = "drop")
site_labels <- with(samples_per_site, setNames(sprintf("%s\n%d / %d", site, n_samples, n_dates), site))

## ---- tests ----
set.seed(42)
tab <- table(cells$site, cells$troph)
tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
overall <- tibble(test = "fisher site x group (cells independent)", group = "all", n = sum(tab),
                  effect = unname(sqrt(chi2 / (sum(tab) * (min(dim(tab)) - 1)))), effect_type = "Cramer's V",
                  p_value = fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value)

glmm_site_p <- function(d) {
  ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  tryCatch({
    m1 <- suppressMessages(suppressWarnings(glmer(y ~ site + (1 | sample), data = d, family = binomial, control = ctrl)))
    m0 <- suppressMessages(suppressWarnings(glmer(y ~ 1 + (1 | sample), data = d, family = binomial, control = ctrl)))
    anova(m0, m1)$`Pr(>Chisq)`[2]
  }, error = function(e) NA_real_)
}
per_group <- map_dfr(group_levels, function(g) {
  d <- cells %>% mutate(y = as.integer(troph == g))
  tibble(test = "GLMM in_group ~ site + (1 | sample), LRT", group = g, n = nrow(d),
         effect = mean(d$y), effect_type = "overall fraction", p_value = glmm_site_p(d))
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
tests <- bind_rows(overall %>% mutate(p_adj = p_value), per_group)

write_tsv(cells, file.path(opt$table_dir, "chloro_foodvac_cells.tsv"))
write_tsv(site_counts, file.path(opt$table_dir, "chloro_foodvac_site_counts.tsv"))
write_tsv(sample_counts, file.path(opt$table_dir, "chloro_foodvac_sample_counts.tsv"))
write_tsv(tests, file.path(opt$table_dir, "chloro_foodvac_tests.tsv"))
print(tests, width = Inf)

## ---- plots ----
fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.0e", p), formatC(p, digits = 2, format = "g"))
glmm_txt <- paste0("GLMM site effect (BH): ",
                   paste0(c("chloro + FV", "chloro only", "no chloro"), " p = ", fmt_p(per_group$p_adj), collapse = ", "))
fisher_txt <- sprintf("Fisher site x group p = %s, Cramer's V = %.2f (cells treated as independent)",
                      fmt_p(overall$p_value), overall$effect)

## segment labels only where the segment is tall enough to hold them
seg_label <- function(n, frac) ifelse(n > 0 & frac >= 0.06, n, "")

p_site <- ggplot(site_counts, aes(site, frac, fill = troph)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = seg_label(n, frac)), position = position_stack(vjust = 0.5),
            colour = "white", size = 3.2) +
  geom_text(data = distinct(site_counts, site, n_site), inherit.aes = FALSE,
            aes(site, 1, label = paste0("n=", n_site)), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = group_cols, name = NULL, drop = FALSE) +
  scale_x_discrete(labels = site_labels) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(x = "site (samples / sampling dates)", y = "cells",
       title = "A  Chloroplast / food vacuole groups per site", subtitle = paste(fisher_txt, glmm_txt, sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), plot.title.position = "plot", legend.position = "top",
        plot.subtitle = element_text(size = 8.5, colour = "grey30"))

sample_counts <- sample_counts %>%
  mutate(sample_label = sprintf("%s %s\n%s", sampling_date, sampling_time, size_frac)) %>%
  arrange(site, sampling_date, sampling_time, sample) %>%
  mutate(sample = factor(sample, unique(sample)))   # samples in date order within each site
p_sample <- ggplot(sample_counts, aes(sample, frac, fill = troph)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = seg_label(n, frac)), position = position_stack(vjust = 0.5),
            colour = "white", size = 2.6) +
  geom_text(data = distinct(sample_counts, site, sample, n_sample), inherit.aes = FALSE,
            aes(sample, 1, label = n_sample), vjust = -0.5, size = 2.6) +
  facet_grid(. ~ site, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_cols, guide = "none", drop = FALSE) +
  scale_x_discrete(labels = setNames(sample_counts$sample_label, sample_counts$sample)) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  labs(x = "sample (sampling date, time, size fraction)", y = "cells",
       title = "B  The same per sample (number on top = cells)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), plot.title.position = "plot",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        strip.text = element_text(face = "bold"))

n_samples <- n_distinct(sample_counts$sample)
save_plot(p_site / p_sample + plot_layout(heights = c(1, 1)), "chloro_foodvac_sites",
          max(8, 0.45 * n_samples + 2), 10, opt$plot_dir)
