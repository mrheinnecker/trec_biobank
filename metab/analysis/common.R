# Shared settings for the metabarcoding analysis scripts (01_ ... 04_).
# Source this at the top of each script. Only this file contains paths.

suppressPackageStartupMessages(library(tidyverse))

# ---- paths -------------------------------------------------------------------
# Works from RStudio on Windows and from WSL/Linux.
base_dir <- if (.Platform$OS.type == "windows") "C:/projects/trec_biobank/metab" else "/mnt/c/projects/trec_biobank/metab"

# Output of the taxonomy workflow (asv_table_taxonomy.tsv, asv_taxonomy.tsv).
# Point this to the new run once it is copied from the cluster.
taxonomy_dir <- file.path(base_dir, "ampliseq_2026-11-03/results/taxonomy")
mapping_f    <- file.path(base_dir, "sample_mapping.tsv")     # from metab/make_sample_mapping.R
# Tara Oceans W4 barcode database (V9, with trophic modes per taxon), for 05_trophic_modes.R
w4_db_f      <- file.path(if (.Platform$OS.type == "windows") "C:/projects" else "/mnt/c/projects",
                          "trec_biobank/analysis/databases/Database_W4_barcode_occurences.tsv")

out_dir   <- file.path(base_dir, "analysis")
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
