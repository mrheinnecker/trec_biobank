## Shared helpers for the TEM annotation analysis.
## Sourced by every numbered script in this folder; not meant to be run on its own.

suppressPackageStartupMessages({
  library(tidyverse)
})

## data locations of all analyses (analysis/paths.R); script_dir is set by the calling script
source(file.path(script_dir, "..", "paths.R"))

## ---- command line ----
## every script takes the same arguments:
##   --outdir      output root; plots go to <outdir>/plots, tables to <outdir>/tables
##   --input       optional: snapshot of the main annotation sheet as tsv
##                 (written by 00_download_annotations.R); without it the main sheet
##                 is read directly from Google Sheets
##   --sheet_url   annotation spreadsheet, --sheet: its tab (default: main)
##   --google_key  service-account key for Google Sheets
##   --rounds      annotation rounds to use, comma separated (default 1 = first annotation only)
## precedence: command line > env vars (TEM_ANNOT_INPUT, TEM_ANNOT_OUTDIR, TEM_ANNOT_SHEET_URL,
##             TEM_ANNOT_SHEET, TEM_GOOGLE_KEY, TEM_ANNOT_ROUNDS) > defaults
## (the defaults are the hard-coded paths at the top of every script, for running line by line)
default_sheet_url  <- "https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=0#gid=0"
default_google_key <- trec_paths$google_key

parse_args <- function(args = commandArgs(trailingOnly = TRUE),
                       default_input = NA_character_, default_outdir = NA_character_) {
  out <- list(input      = Sys.getenv("TEM_ANNOT_INPUT", default_input),
              outdir     = Sys.getenv("TEM_ANNOT_OUTDIR", default_outdir),
              sheet_url  = Sys.getenv("TEM_ANNOT_SHEET_URL", default_sheet_url),
              sheet      = Sys.getenv("TEM_ANNOT_SHEET", "main"),
              google_key = Sys.getenv("TEM_GOOGLE_KEY", default_google_key),
              rounds     = Sys.getenv("TEM_ANNOT_ROUNDS", "1"))
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (str_detect(arg, "^--[^=]+=")) {
      out[[str_remove(str_extract(arg, "^--[^=]+"), "^--")]] <- str_remove(arg, "^--[^=]+=")
      i <- i + 1L
    } else if (str_detect(arg, "^--") && i < length(args)) {
      out[[str_remove(arg, "^--")]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      stop("Unknown or incomplete argument: ", arg)
    }
  }
  if (is.na(out$outdir) || out$outdir == "") {
    stop("Usage: Rscript <script>.R --outdir <output dir> [--input <annotations.tsv>] [--rounds 1]")
  }
  if (!is.na(out$input) && out$input == "") out$input <- NA_character_
  out$rounds <- as.integer(str_split(out$rounds, ",")[[1]])
  if (any(is.na(out$rounds)) || !all(out$rounds %in% 1:2)) stop("--rounds must be 1, 2 or 1,2")
  out$plot_dir  <- file.path(out$outdir, "plots")
  out$table_dir <- file.path(out$outdir, "tables")
  dir.create(out$plot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out$table_dir, showWarnings = FALSE, recursive = TRUE)
  out
}

## ---- output ----
## every plot as pdf + png (white background, theme_minimal is transparent otherwise)
save_plot <- function(p, name, width, height, plot_dir) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(plot_dir, paste0(name, ".", ext)), p, width = width, height = height,
           dpi = 300, bg = "white", limitsize = FALSE)
  }
}

## ---- lookups / feature groups ----
group_lookup <- tibble(
  code  = as.character(1:8),
  group = c("Dinoflagellates", "Diatoms", "Ciliates", "Flagellates",
            "Coccolithophores", "Dictyocha", "Others", "Haptophytes")
)

chloro_feats <- c("chloropl", "plastoglob", "pyrenoid")
common_feats <- c("mito", "ER", "golgi", "nucleus", "nucleolus", "lysosome", "lipid_drop")
meta_feats   <- c("site", "size_frac", "sampling_time", "life_status", "major_group", "beauty")

## organelle / structure columns of the annotation sheet (0/1)
## (listed by name: the column order of the sheet changed with the round format)
organelle_feats <- c(
  "cell_cover", "nucleus", "nucleolus", "tricho", "core_vesic", "flagell_app", "symbiosis",
  "chloropl", "plastoglob", "pyrenoid", "pusule", "retic_net", "electr_sheets", "tubul_net",
  "food_vac", "starch", "cryst_rich_RN", "eyespot", "rhabdo", "ER", "mito", "golgi",
  "lipid_drop", "fibrous_body", "put_division", "SDV", "lysosome", "SER_whorls"
)
organelle_cols <- function(annot) intersect(organelle_feats, names(annot))

## ---- data loading ----
## main annotation sheet (one row per finished annotation round), all columns as text
read_main_table <- function(opt) {
  if (!is.na(opt$input)) {
    message("Annotations from file: ", opt$input)
    return(read_tsv(opt$input, col_types = cols(.default = col_character()), progress = FALSE))
  }
  suppressPackageStartupMessages(library(googlesheets4))
  if (!file.exists(opt$google_key)) {
    stop("Google key not found: ", opt$google_key, " (set --google_key or TEM_GOOGLE_KEY)")
  }
  gs4_auth(path = opt$google_key)
  message("Annotations from Google Sheets, tab '", opt$sheet, "': ", opt$sheet_url)
  read_sheet(opt$sheet_url, sheet = opt$sheet, col_types = "c")
}

## uncertain entries ("?") are counted as absent (0):
##   0/1 structure columns: "??", "1?", "?" -> 0
##   undescribed: numbers marked with "?" ("18?") are dropped, the certain ones stay
##                ("16, 18?" -> "16, "); nothing certain left ("17?", "??") -> "0"
## categorical columns (major_group, taxo_class) are not changed
## every change is kept as attribute "uncertain" of the result (table + figure in 01_annotator_qc.R)
binary_feats <- c(organelle_feats, "thecal_struc", "peduncle")
clean_uncertain <- function(annot) {
  log <- list()
  for (col in intersect(c(binary_feats, "undescribed"), names(annot))) {
    v <- annot[[col]]
    uncertain <- !is.na(v) & str_detect(v, fixed("?"))
    if (!any(uncertain)) next
    cleaned <- if (col == "undescribed") {
      x <- v[uncertain] %>%
        str_remove_all("\\d+\\s*\\?") %>%
        str_remove_all(fixed("?"))
      x <- x %>% str_squish() %>% str_remove_all("^[,;/ ]+|[,;/ ]+$")
      if_else(str_detect(x, "\\d"), x, "0")
    } else {
      rep("0", sum(uncertain))
    }
    log[[col]] <- annot[uncertain, ] %>%
      select(any_of(c("name", "site", "round", "annotated_by"))) %>%
      mutate(column = col, original = v[uncertain], used = cleaned)
    annot[[col]][uncertain] <- cleaned
  }
  log <- bind_rows(log)
  if (nrow(log) > 0) {
    message(sprintf("uncertain entries set to absent: %s",
                    paste0(names(table(log$column)), " ", table(log$column), collapse = ", ")))
  }
  attr(annot, "uncertain") <- log
  annot
}

## annotation rows of the selected rounds (opt$rounds), column types guessed like read_tsv
load_annotations <- function(opt) {
  raw <- read_main_table(opt)
  if (!all(c("name", "round", "annotated_by") %in% names(raw))) {
    stop("The annotation table needs the columns name, round and annotated_by (round format of the main sheet)")
  }
  annot <- raw %>%
    filter(!is.na(name), as.integer(round) %in% opt$rounds) %>%
    clean_uncertain()
  uncertain <- attr(annot, "uncertain")
  annot <- suppressMessages(type_convert(annot))
  attr(annot, "uncertain") <- uncertain
  message(sprintf("%d annotation rows of round %s (%d cells)",
                  nrow(annot), paste(opt$rounds, collapse = ","), n_distinct(annot$name)))
  annot
}

## the first annotation (round 1) of every cell
first_annotation <- function(annot) {
  annot %>% filter(round == 1)
}

## undescribed: messy comma separated numbers ("16, 18?", "9,17,", "7or8", "many 7, 18")
## -> extract every number, one 0/1 column per number (undescribed_<n>)
## "0" means nothing undescribed, NA stays NA in all columns
## uncertain numbers ("17?") were already removed by clean_uncertain() (counted as absent)
add_undescribed_cols <- function(annot) {
  undescribed_wide <- annot %>%
    select(name, undescribed) %>%
    mutate(nr = str_extract_all(undescribed, "\\d+")) %>%
    unnest(nr, keep_empty = TRUE) %>%
    mutate(nr = as.integer(nr)) %>%
    filter(is.na(nr) | nr != 0) %>%
    distinct(name, nr) %>%
    mutate(present = 1L) %>%
    pivot_wider(names_from = nr, values_from = present, values_fill = 0L,
                names_prefix = "undescribed_", names_sort = TRUE) %>%
    select(-any_of("undescribed_NA"))

  annot %>%
    left_join(undescribed_wide, by = "name") %>%
    mutate(across(starts_with("undescribed_"), ~ if_else(is.na(undescribed), NA_integer_, replace_na(.x, 0L))))
}

undescribed_cols <- function(annot) str_subset(names(annot), "^undescribed_\\d+$")

## dinoflagellates (major_group 1), first annotation per cell, with undescribed_<n> columns
## = the data set of the site-feature and association analyses
load_dino_data <- function(opt) {
  load_annotations(opt) %>%
    first_annotation() %>%
    filter(major_group == 1) %>%
    add_undescribed_cols() %>%
    ## I (marco) annotated ret_net also when there is crystaline ret_net ... i am fixing this here
    ## -> retic_net set to 0 wherever cryst_rich_RN is 1
    mutate(
      retic_net = case_when(
        cryst_rich_RN == 1 ~ 0,
        TRUE ~ retic_net
      )
    )
}
