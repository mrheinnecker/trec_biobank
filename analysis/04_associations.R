## Pairwise association tests between all features + metadata (dinoflagellates only)
## test by column type, NAs dropped per pair:
##   0/1 vs 0/1           fisher 2x2, effect = phi, + conditional MLE odds ratio with 95% CI
##   numeric vs numeric   spearman (beauty is ordinal)
##   any categorical      fisher RxC with simulated p, effect = Cramer's V
## Bonferroni over all pairs.
## Output: tables/assoc_results.tsv (+ .rds, read by 05_association_plots.R)

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
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)

cor_data <- build_cor_data(load_dino_data(opt$input))

test_pair <- function(var1, var2, data) {
  d <- data %>% select(x = all_of(var1), y = all_of(var2)) %>% drop_na()
  out <- tibble(var1, var2, n = nrow(d),
                n1 = nrow(d %>% filter(x == 1)),
                n2 = nrow(d %>% filter(y == 1)),
                test = NA_character_, estimate = NA_real_,
                var1_present = NA_character_, var1_absent = NA_character_,
                odds_ratio = NA_real_, or_ci_low = NA_real_, or_ci_high = NA_real_, p_value = NA_real_)
  if (n_distinct(d$x) < 2 || n_distinct(d$y) < 2) return(out)

  if (is_binary(d$x) && is_binary(d$y)) {
    ## 0/1 vs 0/1: exact test on the 2x2 table, effect = phi coefficient (-1..1)
    ## + conditional MLE odds ratio with 95% CI (only defined for 2x2)
    ## counts behind the OR as "var2 present / var2 absent", split by var1 present / absent
    ## OR ~ (a/b) / (c/d)
    tab <- table(factor(d$x, levels = 0:1), factor(d$y, levels = 0:1))
    res <- fisher.test(tab)
    out %>% mutate(test = "fisher 2x2", estimate = cor(d$x, d$y),
                   var1_present = paste0(tab["1", "1"], "/", tab["0", "0"]),
                   var1_absent  = paste0(tab["0", "1"], "/", tab["1", "0"]),
                   x1y1 = tab["1", "1"], x0y0 = tab["0", "0"], x0y1 = tab["0", "1"], x1y0 = tab["1", "0"],
                   odds_ratio = unname(res$estimate), or_ci_low = res$conf.int[1], or_ci_high = res$conf.int[2],
                   p_value = res$p.value)
  } else if (is.numeric(d$x) && is.numeric(d$y)) {
    ## ordinal (beauty) vs numeric: rank correlation (-1..1)
    res <- cor.test(d$x, d$y, method = "spearman", exact = FALSE)
    out %>% mutate(test = "spearman", estimate = unname(res$estimate), p_value = res$p.value)
  } else {
    ## at least one categorical: fisher with simulated p, effect = Cramer's V (0..1)
    tab <- table(d$x, d$y)
    res <- fisher.test(tab, simulate.p.value = TRUE, B = 1e4)
    chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
    out %>% mutate(test = "fisher RxC", estimate = unname(sqrt(chi2 / (sum(tab) * (min(dim(tab)) - 1)))),
                   p_value = res$p.value)
  }
}

set.seed(42)
cor_results <- combn(names(cor_data), 2, simplify = FALSE) %>%
  map_dfr(~ test_pair(.x[1], .x[2], cor_data)) %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni")) %>%
  arrange(p_adj)

write_tsv(cor_results, file.path(opt$table_dir, "assoc_results.tsv"))
#saveRDS(cor_results, file.path(opt$table_dir, "assoc_results.rds"))
message(sprintf("%d pairs tested, %d significant (Bonferroni p < 0.05)",
                sum(!is.na(cor_results$p_value)), sum(cor_results$p_adj < 0.05, na.rm = TRUE)))

