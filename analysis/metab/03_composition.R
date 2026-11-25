# 03: Community composition per sample (QC-passing environmental samples only)
#
# 16S and 18S come from one mixed amplicon, so their read numbers are not
# comparable with each other. Relative abundances are therefore computed
# within each category separately:
#   eukaryote_nuclear  (18S, PR2)        -> share of 18S reads
#   prokaryote         (16S, SILVA)      -> share of prokaryote reads
#   plastid            (16S, PR2 host)   -> share of plastid reads
#
# Output: one stacked-bar plot per category and level (coarse = division /
# phylum, fine = class), plus a site x group table of mean shares.

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/trec_biobank/analysis/metab"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))

asv_info    <- read_prepared("asv_info")
counts_long <- read_prepared("counts_long")
samples     <- read_prepared("samples") %>% filter(pass_qc)

counts <- counts_long %>%
  filter(run_accession %in% samples$run_accession) %>%
  left_join(asv_info %>% select(ASV_ID, category, group_coarse, group_fine), by = "ASV_ID")

sets <- tribble(
  ~category,           ~title,
  "eukaryote_nuclear", "Eukaryotes (18S, PR2)",
  "prokaryote",        "Prokaryotes (16S, SILVA)",
  "plastid",           "Plastids (16S, host lineage from PR2)"
)

composition <- function(cat, level) {
  counts %>%
    filter(category == cat) %>%
    group_by(run_accession, group = .data[[level]]) %>% summarise(reads = sum(reads), .groups = "drop") %>%
    mutate(group = lump_top(group, reads)) %>%
    group_by(run_accession, group) %>% summarise(reads = sum(reads), .groups = "drop") %>%
    group_by(run_accession) %>% mutate(share = reads / sum(reads)) %>% ungroup() %>%
    left_join(samples %>% select(run_accession, site, date, time_of_day, size_fraction, replicate), by = "run_accession")
}

plot_composition <- function(df, title) {
  # Order groups by overall abundance, "other" last.
  lv <- df %>% count(group, wt = reads, sort = TRUE) %>% pull(group)
  df <- df %>% mutate(group = factor(group, levels = c(setdiff(lv, "other"), intersect("other", lv))),
                      sample = paste(date, time_of_day, replicate, run_accession),
                      panel = paste(site, size_fraction, sep = "  "))
  n <- nlevels(df$group)
  pal <- setNames(c(scales::hue_pal(l = 60)(n - ("other" %in% lv)), if ("other" %in% lv) "grey80"), levels(df$group))
  ggplot(df, aes(x = sample, y = share, fill = group)) +
    geom_col(width = 0.9) +
    facet_grid(~ panel, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(x = "samples (ordered by date)", y = "share of reads in this category", fill = NULL, title = title) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0), panel.spacing.x = unit(1, "pt"),
          legend.position = "right")
}

for (i in seq_len(nrow(sets))) {
  for (level in c("group_coarse", "group_fine")) {
    df <- composition(sets$category[i], level)
    if (nrow(df) == 0) next
    lvl_name <- if (level == "group_coarse") "coarse" else "fine"
    save_plot(plot_composition(df, paste0(sets$title[i], " - ", lvl_name)),
              paste0("03_composition_", sets$category[i], "_", lvl_name), width = 16, height = 6)

    # Mean share per site and size fraction (samples weighted equally).
    df %>% group_by(site, size_fraction, group) %>%
      summarise(mean_share = sum(share) / n_distinct(run_accession), .groups = "drop") %>%
      complete(nesting(site, size_fraction), group, fill = list(mean_share = 0)) %>%
      pivot_wider(names_from = group, values_from = mean_share) %>%
      write_tsv(file.path(table_dir, paste0("03_mean_share_", sets$category[i], "_", lvl_name, ".tsv")))
  }
}

# Plastid host lineage needs the organelle step of the taxonomy workflow.
if (all(is.na(asv_info$organelle_host_taxonomy)))
  message("No organelle host lineage in this taxonomy run: plastids are all 'host unknown'. ",
          "Rerun with the new taxonomy output to get them.")

cat("plots written to", plot_dir, "\n")
