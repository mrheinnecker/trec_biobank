## Organelles x organelles: pairwise association tests between the selected organelles
## (use_features + undescribed_<n>, set in assoc_common.R; dinoflagellates only)
## Only pairs of two organelles; the organelles vs metadata are a separate analysis (07/08).
## Tests (see test_pair() in assoc_common.R): fisher 2x2, effect = phi + odds ratio with 95% CI
## Bonferroni over the pairs of this analysis.
## Output: tables/assoc_organelles_results.tsv (+ .rds, read by 05_organelle_association_plots.R)

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

dino <- load_dino_data(opt)
cor_data <- build_cor_data(dino)
feats <- assoc_features(dino, verbose = TRUE)

pairs <- combn(feats, 2, simplify = FALSE) %>%
  map_dfr(~ tibble(var1 = .x[1], var2 = .x[2]))

results <- run_assoc_tests(pairs, cor_data)
save_assoc_results(results, "assoc_organelles", opt$table_dir)
