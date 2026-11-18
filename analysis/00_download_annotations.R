## Download the main annotation sheet from Google Sheets into one snapshot,
## so that all analysis steps of a run work on the same data.
## Writes <outdir>/tables/annotations_main.tsv (used by run_analysis.sh as --input for the
## other steps) and a timestamped copy annotations_main_<YYYYmmdd_HHMMSS>.tsv.
## All rounds are kept in the snapshot; the steps select the rounds with --rounds.

## ---- paths ----
default_outdir     <- "C:/projects/trec_biobank"
default_script_dir <- "C:/repos/tem_classification/analysis"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
opt <- parse_args(default_outdir = default_outdir)

if (!is.na(opt$input)) stop("00_download_annotations.R reads from Google Sheets; do not pass --input")

main <- read_main_table(opt)
if (!all(c("name", "round", "annotated_by") %in% names(main))) {
  stop("The main sheet has no name/round/annotated_by columns (not in the round format)")
}

snapshot <- file.path(opt$table_dir, "annotations_main.tsv")
write_tsv(main, snapshot, na = "")
invisible(file.copy(snapshot, file.path(opt$table_dir, paste0("annotations_main_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"))))
message(sprintf("%d rows (%s) -> %s", nrow(main),
                paste0("round ", names(table(main$round)), ": ", table(main$round), collapse = ", "),
                snapshot))
