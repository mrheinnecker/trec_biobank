## Shared helpers for the TEM annotation analysis.
## Sourced by every numbered script in this folder; not meant to be run on its own.

suppressPackageStartupMessages({
  library(tidyverse)
})

## ---- command line ----
## every script takes the same arguments:
##   --input   annotation table (tsv, export of the annotation sheet)
##   --outdir  output root; plots go to <outdir>/plots, tables to <outdir>/tables
## precedence: command line > env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR > defaults
## (the defaults are the hard-coded paths at the top of every script, for running line by line)
parse_args <- function(args = commandArgs(trailingOnly = TRUE),
                       default_input = NA_character_, default_outdir = NA_character_) {
  out <- list(input  = Sys.getenv("TEM_ANNOT_INPUT", default_input),
              outdir = Sys.getenv("TEM_ANNOT_OUTDIR", default_outdir))
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
  if (is.na(out$input) || is.na(out$outdir)) {
    stop("Usage: Rscript <script>.R --input <annotations.tsv> --outdir <output dir>")
  }
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
organelle_cols <- function(annot) names(select(annot, cell_cover:SER_whorls))

## ---- data loading ----
load_annotations <- function(path) {
  read_tsv(path, show_col_types = FALSE)
}

## cells annotated in several runs: keep only the annotation from the earliest run
## run columns are YYYYMMDD_HHMMSS timestamps -> chronological = alphabetical
dedup_first_run <- function(annot) {
  first_run <- annot %>%
    filter(!is.na(file)) %>%
    group_by(name) %>%
    summarise(run = min(run), n_runs = n(), .groups = "drop")
  annot %>% semi_join(first_run, by = c("name", "run"))
}

## undescribed: messy comma separated numbers ("16, 18?", "9,17,", "7or8", "many 7, 18")
## -> extract every number, one 0/1 column per number (undescribed_<n>)
## "0" means nothing undescribed, NA stays NA in all columns
## uncertain entries ("17?") are counted as present
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
load_dino_data <- function(path) {
  load_annotations(path) %>%
    dedup_first_run() %>%
    filter(major_group == 1) %>%
    add_undescribed_cols() %>%
    ## I (marco) annotated ret_net also when there is crystaline ret_net ... i am fixing this here
    mutate(
      retic_net_new=case_when(
        cryst_rich_RN==1 ~ 0,
        TRUE ~ retic_net
      )
    ) #%>% select(retic_net, cryst_rich_RN,  retic_net_new) %>% View()
}
