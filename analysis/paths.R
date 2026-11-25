## Data locations of all analyses (tem/, metab/, overview/). Only this file holds the paths.
## Sourced by tem/common.R, metab/common.R and the overview scripts.
## Every entry can be overridden with an env var (set by analysis/launcher.sh); the defaults are for
## running a script by hand in RStudio (Windows) or with Rscript in WSL.
##
## <data_root>/                       TREC_DATA_ROOT
##   analysis/
##     databases/                     shared reference tables (Tara W4 barcode database)
##     tem/<run>/                     TEM annotation analysis, one folder per run
##     metab/<run>/                   metabarcoding analysis, one folder per run
##     overview/data/<snapshot>/      biobank overview: data snapshots
##     overview/runs/<run>/           biobank overview: figures
##     overview/event_overrides.tsv   manual event mappings
##   metab/                           metabarcoding pipeline output (ampliseq, taxonomy, sample_mapping.tsv)
##   pngs/, png_thumbs/               TEM overview PNGs (wfTEM) and their reduced copies

trec_env <- function(var, default) { v <- Sys.getenv(var); if (nzchar(v)) v else default }

trec_data_root <- trec_env("TREC_DATA_ROOT",
                           if (.Platform$OS.type == "windows") "C:/projects/trec_biobank" else "/mnt/c/projects/trec_biobank")
trec_repo      <- trec_env("TREC_REPO",
                           if (.Platform$OS.type == "windows") "C:/repos/trec_biobank" else "/mnt/c/repos/trec_biobank")

trec_paths <- local({
  analysis <- file.path(trec_data_root, "analysis")
  databases <- file.path(analysis, "databases")
  list(
    data_root      = trec_data_root,
    analysis       = analysis,
    databases      = databases,
    ## outputs, one folder per run inside
    tem_runs       = trec_env("TREC_TEM_RUNS", file.path(analysis, "tem")),
    metab_runs     = trec_env("TREC_METAB_RUNS", file.path(analysis, "metab")),
    overview       = trec_env("TREC_OVERVIEW_ROOT", file.path(analysis, "overview")),
    ## inputs
    w4_db          = trec_env("METAB_W4_DB", file.path(databases, "Database_W4_barcode_occurences.tsv")),
    metab_dir      = file.path(trec_data_root, "metab"),
    metab_taxonomy = trec_env("METAB_TAXONOMY_DIR", file.path(trec_data_root, "metab/ampliseq_2026-11-19/taxonomy")),
    metab_mapping  = trec_env("METAB_MAPPING", file.path(trec_data_root, "metab/sample_mapping.tsv")),
    png_dir        = trec_env("TEM_PNG_DIR", file.path(trec_data_root, "pngs")),
    thumb_dir      = trec_env("TEM_THUMB_DIR", file.path(trec_data_root, "png_thumbs")),
    google_key     = trec_env("TREC_GOOGLE_KEY", file.path(trec_repo, "trec-tem-screen-e98a2e03f58b.json"))
  )
})

## newest run folder of a module that contains `file` (e.g. the metab table that tem/13 needs)
trec_latest_run <- function(runs_dir, file = NULL) {
  runs <- sort(list.dirs(runs_dir, recursive = FALSE), decreasing = TRUE)
  if (!is.null(file)) runs <- runs[file.exists(file.path(runs, file))]
  if (length(runs) == 0) NA_character_ else runs[1]
}
