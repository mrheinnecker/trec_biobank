## Organelles x metadata: every selected organelle (use_features + undescribed_<n>, set in
## assoc_common.R) tested against the metadata (site, size_frac, sampling_time, life_status,
## beauty; dinoflagellates only). Separate from the organelle x organelle analysis (04/05),
## with its own Bonferroni correction.
## Tests (see test_pair() in assoc_common.R):
##   organelle vs site / size_frac / sampling_time / life_status: fisher RxC (simulated p), Cramer's V
##   organelle vs beauty (ordinal): spearman
## Output: tables/assoc_metadata_results.tsv (+ .rds, read by 08_metadata_association_plots.R)

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
feats <- assoc_features(dino)
meta <- assoc_metadata(cor_data)
message("metadata: ", paste(meta, collapse = ", "))

## var1 = organelle, var2 = metadata
pairs <- expand_grid(var1 = feats, var2 = meta)

results <- run_assoc_tests(pairs, cor_data)
save_assoc_results(results, "assoc_metadata", opt$table_dir)
