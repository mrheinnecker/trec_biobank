# Shared settings for the metabarcoding analysis scripts (01_ ... 04_).
# Source this at the top of each script. Only this file contains paths.

suppressPackageStartupMessages(library(tidyverse))

# ---- paths -------------------------------------------------------------------
# All data locations: analysis/paths.R (script_dir is set by the calling script).
# The defaults there are for running a script by hand in RStudio; analysis/launcher.sh sets the
# METAB_* variables instead, with one output folder per run.
source(file.path(script_dir, "..", "paths.R"))
env_or <- trec_env

taxonomy_dir <- trec_paths$metab_taxonomy   # output of the taxonomy workflow (asv_table_taxonomy.tsv, asv_taxonomy.tsv)
mapping_f    <- trec_paths$metab_mapping    # from metab/make_sample_mapping.R
w4_db_f      <- trec_paths$w4_db            # Tara Oceans W4 barcode database (V9, trophic modes), for 05_trophic_modes.R

out_dir   <- env_or("METAB_OUTDIR", file.path(trec_paths$metab_runs, "rstudio"))
data_dir  <- file.path(out_dir, "data")      # prepared R objects (01_)
table_dir <- file.path(out_dir, "tables")
plot_dir  <- file.path(out_dir, "plots")
for (d in c(data_dir, table_dir, plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- analysis settings -------------------------------------------------------
min_sample_reads <- 10000   # samples with fewer ASV reads are excluded downstream
top_n_groups     <- 12      # taxa shown individually in composition plots, rest = "other"

# ---- helpers -----------------------------------------------------------------
read_prepared <- function(name) readRDS(file.path(data_dir, paste0(name, ".rds")))

save_plot <- function(p, name, width = 10, height = 6) {
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(plot_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 150)
}

# Collapse all but the n most abundant levels of `group` (by total `reads`) into "other".
lump_top <- function(group, reads, n = top_n_groups) {
  top <- tibble(group, reads) %>% filter(!is.na(group)) %>% count(group, wt = reads, sort = TRUE) %>%
    slice_head(n = n) %>% pull(group)
  if_else(group %in% top, group, "other")
}

category_colours <- c(prokaryote = "#4C72B0", plastid = "#55A868", mitochondrion = "#C44E52",
                      eukaryote_nuclear = "#DD8452", unassigned = "grey70")

theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
                                           strip.background = element_rect(fill = "grey95")))
