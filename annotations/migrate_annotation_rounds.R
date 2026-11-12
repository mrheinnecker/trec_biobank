## One-off migration of the TEM annotations into the round format.
##
## Rebuilds the main table (one row per finished round) and the status table
## (one row per image: annotated_by, validated_by) from all backups in the
## annotation log dir, because in the old format the validation overwrote the
## first annotation in the main sheet.
##
## Sources, per backup folder <log_dir>/<timestamp>/:
##   <person>.tsv   personal tabs (main source)
##   main sheet     (annotations_sheet from run_configuration.tsv) only as fallback for
##                  image/person combinations that are in no personal tab, and only rows
##                  without validated_by (validated rows of the old main are mixed)
## In google mode the live sheets are first backed up into a new folder, so work
## done since the last run is included.
##
## Rules (see resolve_rounds() in annotations/annotation_rounds_common.R):
##   - the name in annotated_by / validated_by counts, not the tab it is in
##   - rows with at least one annotation field filled and a name count as finished
##   - one version per image and person: most filled fields, then the newest backup
##   - validated_by -> round 2; otherwise the first person (oldest backup) is round 1,
##     the next different person round 2
##
## Without --write only TSVs + an issues report are written to --outdir.
## With --write TRUE the main and status sheets are overwritten (personal tabs are not touched).
##
## Example (cluster, tidyverse container):
## Rscript annotations/migrate_annotation_rounds.R \
##   --google_key /g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json \
##   --outdir /g/schwab/tem_screen/annotations/migration_rounds \
##   --write FALSE

library(tidyverse)
library(getopt)

spec <- matrix(c(
  "log_dir", "a", 1, "character",
  "outdir", "o", 1, "character",
  "assignees", "g", 1, "character",
  "sheet_mode", "m", 1, "character",
  "google_key", "k", 1, "character",
  "collection_table_url", "u", 1, "character",
  "source_collection_table_url", "c", 1, "character",
  "source_collection_sheet", "t", 1, "character",
  "local_source_collection_table", "f", 1, "character",
  "annotations_sheet", "n", 1, "character",
  "status_sheet", "s", 1, "character",
  "write", "w", 1, "character"
), ncol=4, byrow=TRUE)
opt <- getopt(spec)

arg_file <- commandArgs(FALSE) %>%
  .[str_detect(., "^--file=")] %>%
  str_remove("^--file=") %>%
  .[1]
script_dir <- if (!is.na(arg_file)) dirname(normalizePath(arg_file)) else getwd()
source(file.path(script_dir, "annotation_rounds_common.R"))

opt_or <- function(value, default) {
  if (is.null(value) || is.na(value)) default else value
}

log_dir <- opt_or(opt$log_dir, "/g/schwab/tem_screen/annotations/log")
outdir <- opt_or(opt$outdir, "migration_rounds")
people <- str_split(opt_or(opt$assignees, "marco,chandni,yannick,karel,viktoria"), ",")[[1]] %>%
  str_trim() %>%
  str_to_lower() %>%
  discard(~.x == "")
sheet_mode <- opt_or(opt$sheet_mode, "google")
google_key <- opt_or(opt$google_key, "/g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json")
collection_table_url <- opt_or(
  opt$collection_table_url,
  "https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=0#gid=0"
)
source_collection_table_url <- opt_or(
  opt$source_collection_table_url,
  "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=2014112369#gid=2014112369"
)
source_collection_sheet <- opt_or(opt$source_collection_sheet, "tem_collection_table")
local_source_collection_table <- opt_or(opt$local_source_collection_table, "tem_collection_table.tsv")
annotations_sheet <- opt_or(opt$annotations_sheet, "main")
status_sheet <- opt_or(opt$status_sheet, "status")
do_write <- toupper(opt_or(opt$write, "FALSE")) %in% c("TRUE", "T", "1", "YES")

if (do_write && sheet_mode != "google") {
  stop("--write TRUE only works with --sheet_mode google")
}
if (!dir.exists(log_dir)) {
  stop("Log dir does not exist: ", log_dir)
}
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

read_chr_tsv <- function(path) {
  if (file.size(path) == 0) {
    return(tibble())
  }
  read_tsv(path, col_types=cols(.default=col_character()), progress=FALSE)
}

## ---- 1. back up the live sheets (google mode) ----
if (sheet_mode == "google") {
  library(googlesheets4)
  library(googledrive)
  gs4_auth(path=google_key)
  drive_auth(path=google_key)

  live_dir <- file.path(log_dir, format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (dir.exists(live_dir)) stop("Backup directory already exists: ", live_dir)
  dir.create(live_dir, recursive=TRUE)
  existing_sheets <- sheet_properties(collection_table_url)$name
  for (sheet in c(annotations_sheet, people)) {
    table <- if (sheet %in% existing_sheets) {
      read_sheet(collection_table_url, sheet=sheet, col_types="c")
    } else {
      tibble()
    }
    if (ncol(table) == 0) write_lines(character(), file.path(live_dir, paste0(sheet, ".tsv")))
    else write_tsv(table, file.path(live_dir, paste0(sheet, ".tsv")), na="")
  }
  write_tsv(
    tibble(
      setting=c("sheet_mode", "collection_table_url", "annotations_sheet", "assignees", "max_rows_per_person"),
      value=c(sheet_mode, collection_table_url, annotations_sheet, paste(people, collapse=","), "migration")
    ),
    file.path(live_dir, "run_configuration.tsv")
  )
  message("Live sheets backed up to: ", live_dir)
}

## ---- 2. read all backups ----
backup_dirs <- list.dirs(log_dir, recursive=FALSE, full.names=TRUE) %>%
  keep(~str_detect(basename(.x), "^\\d{8}_\\d{6}$")) %>%
  sort()
if (length(backup_dirs) == 0) {
  stop("No backup folders (YYYYmmdd_HHMMSS) in ", log_dir)
}

not_tab_files <- c("merged_annotations.tsv", "run_configuration.tsv", "issues.tsv")

backup_files <- map_dfr(backup_dirs, function(dir) {
  config_file <- file.path(dir, "run_configuration.tsv")
  main_name <- "annotations_log"
  if (file.exists(config_file)) {
    config <- read_chr_tsv(config_file)
    value <- config$value[config$setting == "annotations_sheet"]
    if (length(value) == 1 && !is.na(value)) main_name <- value
  }
  files <- list.files(dir, pattern="\\.tsv$", full.names=TRUE)
  tibble(path=files, file=basename(files), run=basename(dir)) %>%
    filter(!file %in% not_tab_files, !str_starts(file, "new_"), file != "status.tsv") %>%
    mutate(kind=if_else(file %in% paste0(c(main_name, "main", "annotations_log"), ".tsv"), "main", "tab"))
})

read_backup <- function(files) {
  pmap_dfr(files, function(path, file, run, kind) {
    table <- read_chr_tsv(path)
    if (nrow(table) == 0 || !"name" %in% names(table)) {
      return(tibble())
    }
    if (kind == "main" && "round" %in% names(table)) {
      message("Skipping main backup in round format: ", path)
      return(tibble())
    }
    table %>%
      mutate(.sheet=str_remove(file, "\\.tsv$"), .run=run, .source=path)
  })
}

tab_rows <- classify_tab_rows(read_backup(backup_files %>% filter(kind == "tab")))
main_rows <- classify_tab_rows(read_backup(backup_files %>% filter(kind == "main")))
message("Backup folders: ", length(backup_dirs), " | tab rows: ", nrow(tab_rows),
        " | old main rows: ", nrow(main_rows))

issues <- list()

## values entered without any name
unnamed <- tab_rows %>% filter(.status == "in_progress")
if (nrow(unnamed) > 0) {
  issues$unnamed <- unnamed %>%
    transmute(type="values_without_name_ignored", name, person=.sheet,
              detail=paste0(.n_filled, " fields filled, no annotated_by/validated_by; ", .source))
}

tab_candidates <- tab_rows %>% filter(.status == "signed")

## fallback from the old main sheet: signed, not validated, not in any tab
main_candidates <- main_rows %>%
  filter(.status == "signed", is.na(.val)) %>%
  anti_join(tab_candidates %>% distinct(name, .person), by=c("name", ".person"))
if (nrow(main_candidates) > 0) {
  issues$main_fallback <- main_candidates %>%
    group_by(name, .person) %>%
    summarise(detail=paste0("only found in old main sheet (", n(), " backups)"), .groups="drop") %>%
    transmute(type="taken_from_old_main_sheet", name, person=.person, detail)
}

candidates <- bind_rows(tab_candidates, main_candidates)

unknown <- candidates %>% filter(!.person %in% people)
if (nrow(unknown) > 0) {
  issues$unknown <- unknown %>%
    distinct(name, .person, .source) %>%
    transmute(type="unknown_annotator_name", name, person=.person, detail=.source)
}

## ---- 3. rebuild the rounds ----
resolved <- resolve_rounds(tibble(), candidates)

## ---- 4. collection table, identity columns, status ----
collection <- if (sheet_mode == "google") {
  read_sheet(source_collection_table_url, sheet=source_collection_sheet, col_types="c")
} else {
  read_chr_tsv(local_source_collection_table)
}
collection <- collection %>%
  filter(!is_blank(name)) %>%
  distinct(name, .keep_all=TRUE) %>%
  select(any_of(collection_identity_columns))

rounds <- resolved$rounds %>%
  attach_identity(collection) %>%
  sort_rounds()

status <- make_status_table(collection, rounds)

not_in_collection <- setdiff(unique(rounds$name), collection$name)
if (length(not_in_collection) > 0) {
  issues$not_in_collection <- tibble(
    type="not_in_collection_table", name=not_in_collection, person=NA_character_,
    detail="annotated image is not in the TEM collection table (kept in main)"
  )
}

low <- rounds %>%
  mutate(.n_filled=count_filled(rounds)) %>%
  filter(.n_filled <= 5)
if (nrow(low) > 0) {
  issues$low <- low %>%
    transmute(type="few_fields_filled", name, person=annotated_by,
              detail=paste0("round ", round, ": ", .n_filled, " fields filled"))
}

issues <- bind_rows(resolved$issues, bind_rows(issues)) %>%
  arrange(type, name)

## ---- 5. output ----
write_tsv(rounds, file.path(outdir, paste0(annotations_sheet, "_rounds.tsv")), na="")
write_tsv(status, file.path(outdir, paste0(status_sheet, ".tsv")), na="")
write_tsv(issues, file.path(outdir, "issues.tsv"), na="")

message("Rounds: ", nrow(rounds), " (round 1: ", sum(rounds$round == 1L),
        ", round 2: ", sum(rounds$round == 2L), ") for ", n_distinct(rounds$name), " images")
message("Collection table images: ", nrow(collection),
        " | with round 1: ", sum(!is.na(status$annotated_by)),
        " | with round 2: ", sum(!is.na(status$validated_by)))
message("Issues:")
print(issues %>% count(type), n=Inf)
message("Written to: ", normalizePath(outdir))

if (do_write) {
  if (any(issues$type == "unknown_annotator_name")) {
    stop("Unknown annotator names (see issues.tsv). Fix them or add them to name_aliases; nothing was written.")
  }
  write_sheet(rounds, ss=collection_table_url, sheet=annotations_sheet)
  write_sheet(status, ss=collection_table_url, sheet=status_sheet)
  message("Sheets '", annotations_sheet, "' and '", status_sheet, "' overwritten.")
}
