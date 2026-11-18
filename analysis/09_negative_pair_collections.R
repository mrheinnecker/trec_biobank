## MoBIE collection tables for the negatively associated organelle pairs (from 04_organelle_associations.R)
## For every significant pair with a negative association (fisher 2x2, p_adj < 0.05, phi < 0) all cells
## of the analysis (dinoflagellates, selected rounds, both organelles annotated) are split into
##   <A>_only, <B>_only, both, none
## and every group gets its own view/grid, so MoBIE opens each group as a separate grid view.
## Output:
##   tables/collections/<A>_vs_<B>.tsv          one collection table per pair
##   tables/collections/overview.tsv            pairs with the number of cells per group
##   Google spreadsheet (collection_sheet_url): one tab per pair + tab "overview"
##   (the service account needs edit access to that spreadsheet; tabs of earlier runs that are no
##    longer a negative pair are not deleted, see the overview tab for the current set)

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

## ---- settings ----
## spreadsheet: TEM_COLLECTION_SHEET_URL (collection_sheet() in assoc_common.R); columns: collection_cols
p_cutoff <- 0.05

## ---- data ----
dino <- load_dino_data(opt)
results <- read_assoc_results("assoc_organelles", opt$table_dir, "04_organelle_associations.R")

neg_pairs <- results %>%
  filter(test == "fisher 2x2", p_adj < p_cutoff, estimate < 0) %>%
  arrange(p_adj)
if (nrow(neg_pairs) == 0) {
  message("No significant negative organelle pairs, no collection tables")
  if (!interactive()) quit(save = "no")
}

## ---- one collection table per pair ----
make_collection <- function(a, b) {
  dino %>%
    filter(!is.na(.data[[a]]), !is.na(.data[[b]])) %>%
    mutate(
      group = case_when(
        .data[[a]] == 1 & .data[[b]] == 0 ~ paste0(a, "_only"),
        .data[[a]] == 0 & .data[[b]] == 1 ~ paste0(b, "_only"),
        .data[[a]] == 1 & .data[[b]] == 1 ~ "both",
        TRUE                              ~ "none"
      ),
      group = factor(group, c(paste0(a, "_only"), paste0(b, "_only"), "both", "none")),
      view = as.character(group),
      grid = as.character(group)
    ) %>%
    arrange(group, site, name) %>%
    select(any_of(collection_cols), group, all_of(c(a, b)), annotated_by, round) %>%
    mutate(group = as.character(group))
}

## tab names: "<A>_vs_<B>" (Google Sheets allows up to 100 characters)
collections <- neg_pairs %>%
  mutate(tab = str_trunc(paste0(var1, "_vs_", var2), 100, ellipsis = "")) %>%
  select(tab, var1, var2, estimate, p_adj) %>%
  mutate(table = map2(var1, var2, make_collection))

n_of <- function(counts, key) { v <- unname(counts[key]); if (is.na(v)) 0L else as.integer(v) }
overview <- collections %>%
  mutate(counts = map(table, ~ count(.x, group) %>% pull(n, group))) %>%
  transmute(
    tab, organelle_a = var1, organelle_b = var2, phi = round(estimate, 3), p_adj = signif(p_adj, 3),
    a_only  = map2_int(counts, organelle_a, ~ n_of(.x, paste0(.y, "_only"))),
    b_only  = map2_int(counts, organelle_b, ~ n_of(.x, paste0(.y, "_only"))),
    both    = map_int(counts, ~ n_of(.x, "both")),
    none    = map_int(counts, ~ n_of(.x, "none"))
  )

collection_dir <- file.path(opt$table_dir, "collections")
dir.create(collection_dir, showWarnings = FALSE, recursive = TRUE)
walk2(collections$table, collections$tab, ~ write_tsv(.x, file.path(collection_dir, paste0(.y, ".tsv")), na = ""))
write_tsv(overview, file.path(collection_dir, "overview.tsv"))
message(sprintf("%d negative pairs -> %s", nrow(collections), collection_dir))
print(overview, n = Inf, width = Inf)

## ---- push to Google Sheets ----
push_collection_tabs(c(list(overview = overview), set_names(collections$table, collections$tab)), opt)
