# 05: Trophic modes of the 18S ASVs from the Tara Oceans W4 database (de Vargas et al. 2015)
#
# The W4 database is 18S V9, our amplicon is 18S V4-V5, so ASVs cannot be matched by
# sequence. W4 assigns one trophic mode per taxon (genus), so the ASVs are matched by
# their PR2 names instead:
#   1. PR2 species (genus_epithet) against the W4 species
#   2. else PR2 genus against the W4 genus
#   W4 genus names carry clade numbers (Gyrodinium_05) -> removed; old names translated
#   via `w4_synonyms` (Neoceratium -> Tripos). No match above genus level (trophic modes
#   differ within families). No manual corrections.
# W4 trophic modes: phototroph (has a plastid, mixotrophs included), heterotroph,
#   photosymbiosis (hosts of photosymbionts / kleptoplastidy), Metazoa.
# ASVs without a mode are split into
#   "NA: genus not in W4"   PR2 genus known, but W4 has no trophic mode for it
#   "NA: no genus (PR2)"    PR2 did not assign a genus (or only a placeholder like Dinophyceae_XXX)
#
# All samples of the count table are used (no QC / site filter); shares within the
# 18S (eukaryote_nuclear) reads of each sample.
# Output:
#   tables/05_asv_trophic_modes.tsv        one row per 18S ASV: taxonomy, reads, match, mode
#   tables/05_trophic_mode_per_sample.tsv  reads and share per sample x mode (all eukaryotes / dinoflagellates)
#   plots/05_trophic_modes_eukaryotes, plots/05_trophic_modes_dinoflagellates

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/tem_classification/metab/analysis"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))

w4_synonyms <- c(neoceratium = "tripos")      # W4 name -> PR2 name (lower case)

mode_levels <- c("phototroph", "photosymbiosis", "heterotroph", "Metazoa",
                 "NA: genus not in W4", "NA: no genus (PR2)")
mode_colours <- c(phototroph = "#1B9E77", photosymbiosis = "#A6D854", heterotroph = "#D95F02",
                  Metazoa = "#5B7DB1", `NA: genus not in W4` = "grey55", `NA: no genus (PR2)` = "grey82")

# ---- 1. W4 lookup: taxon -> trophic mode (cached) ----------------------------
norm_taxon <- function(x) {
  x <- str_remove(str_replace_all(tolower(x), "[+ ]", "_"), "_\\d+$")
  if_else(x %in% names(w4_synonyms), unname(w4_synonyms[x]), x)
}

lookup_f <- file.path(data_dir, "w4_trophic_lookup.rds")
if (!file.exists(lookup_f) || file.mtime(lookup_f) < file.mtime(w4_db_f)) {
  cat("reading", w4_db_f, "(only lineage + trophism, takes a few minutes)\n")
  w4 <- data.table::fread(w4_db_f, select = c("lineage", "trophism"), na.strings = c("", "NA"))
  lin <- as_tibble(w4[!is.na(trophism), .N, by = .(lineage, trophism)])
  rm(w4); invisible(gc())
  # a genus token (capitalised, optional clade number) followed by a lower-case epithet = species
  tok <- lin %>%
    mutate(id = row_number(), t = str_split(lineage, "\\|")) %>%
    unnest(t) %>%
    group_by(id) %>% mutate(nxt = lead(t)) %>% ungroup() %>%
    mutate(is_genus = str_detect(t, "^[A-Z][a-z]+(_\\d+)?$"),
           species  = if_else(is_genus & !is.na(nxt) & str_detect(nxt, "^[a-z]+$") & nxt != "sp",
                              paste0(norm_taxon(t), "_", nxt), NA_character_))
  w4_lookup <- bind_rows(
    tok %>% filter(is_genus)       %>% transmute(rank = "genus",   taxon = norm_taxon(t), id, trophism, N),
    tok %>% filter(!is.na(species)) %>% transmute(rank = "species", taxon = species,        id, trophism, N)) %>%
    distinct() %>%
    group_by(rank, taxon, trophism) %>% summarise(n = sum(N), .groups = "drop_last") %>%
    mutate(w4_barcodes = sum(n), purity = n / w4_barcodes) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(rank, taxon, trophism, purity, w4_barcodes)
  saveRDS(w4_lookup, lookup_f)
} else {
  w4_lookup <- readRDS(lookup_f)
}
cat("W4 lookup:", sum(w4_lookup$rank == "genus"), "genera,", sum(w4_lookup$rank == "species"), "species\n")

# ---- 2. ASVs -> trophic mode -------------------------------------------------
asv_info    <- read_prepared("asv_info")
counts_long <- read_prepared("counts_long")
samples     <- read_prepared("samples") %>% filter(in_count_table)

euk <- asv_info %>% filter(category == "eukaryote_nuclear") %>%
  select(ASV_ID, starts_with("pr2_"))

placeholder <- function(x) is.na(x) | str_detect(x, "_X|XX|_sp\\.?$")
matches <- bind_rows(
  euk %>% filter(!placeholder(pr2_Species)) %>% transmute(ASV_ID, rank = "species", taxon = norm_taxon(pr2_Species)),
  euk %>% filter(!placeholder(pr2_Genus))   %>% transmute(ASV_ID, rank = "genus",   taxon = norm_taxon(pr2_Genus))) %>%
  inner_join(w4_lookup, by = c("rank", "taxon")) %>%
  mutate(rank = factor(rank, c("species", "genus"))) %>%
  group_by(ASV_ID) %>% slice_min(as.integer(rank), n = 1, with_ties = FALSE) %>% ungroup()

asv_reads <- counts_long %>% filter(ASV_ID %in% euk$ASV_ID) %>% count(ASV_ID, wt = reads, name = "reads")
asv_modes <- euk %>%
  left_join(matches %>% select(ASV_ID, match_rank = rank, match_taxon = taxon, trophism, purity, w4_barcodes), by = "ASV_ID") %>%
  left_join(asv_reads, by = "ASV_ID") %>%
  mutate(reads = replace_na(reads, 0),
         mode = case_when(!is.na(trophism)        ~ trophism,
                          !placeholder(pr2_Genus) ~ "NA: genus not in W4",
                          TRUE                    ~ "NA: no genus (PR2)"),
         mode = factor(mode, mode_levels))
write_tsv(asv_modes, file.path(table_dir, "05_asv_trophic_modes.tsv"), na = "")

cat("\nreads per trophic mode, all eukaryotes:\n")
print(asv_modes %>% count(mode, wt = reads, name = "reads") %>% mutate(share = round(reads / sum(reads), 3)))

# ---- 3. Share per sample -----------------------------------------------------
per_sample <- function(asv_set, set_name) {
  counts_long %>%
    filter(run_accession %in% samples$run_accession) %>%
    inner_join(asv_set %>% select(ASV_ID, mode), by = "ASV_ID") %>%
    count(run_accession, mode, wt = reads, name = "reads") %>%
    complete(run_accession, mode, fill = list(reads = 0)) %>%
    group_by(run_accession) %>% mutate(set_reads = sum(reads), share = reads / set_reads) %>% ungroup() %>%
    filter(set_reads > 0) %>%
    mutate(set = set_name)
}
shares <- bind_rows(per_sample(asv_modes, "all eukaryotes"),
                    per_sample(asv_modes %>% filter(pr2_Class == "Dinophyceae"), "dinoflagellates")) %>%
  left_join(samples %>% select(run_accession, sample_type, site, date, time_of_day, size_fraction, replicate, label),
            by = "run_accession")
write_tsv(shares, file.path(table_dir, "05_trophic_mode_per_sample.tsv"), na = "")

# ---- 4. Plots ----------------------------------------------------------------
plot_modes <- function(df, title) {
  df <- df %>%
    mutate(panel  = case_when(sample_type == "control" ~ "controls",
                              is.na(site)              ~ "no metadata",
                              TRUE                     ~ paste(site, coalesce(size_fraction, "?"), sep = "  ")),
           sample = paste(coalesce(as.character(date), ""), coalesce(time_of_day, ""), coalesce(replicate, ""), run_accession)) %>%
    arrange(panel, date, sample) %>%
    mutate(sample = factor(sample, unique(sample)))
  n_na <- df %>% filter(str_starts(mode, "NA")) %>% group_by(run_accession) %>% summarise(s = sum(share)) %>% pull(s)
  ggplot(df, aes(sample, share, fill = mode)) +
    geom_col(width = 0.9) +
    facet_grid(~ panel, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = mode_colours, drop = FALSE) +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.01))) +
    labs(x = "samples (ordered by date)", y = "share of reads", fill = "trophic mode (W4)", title = title,
         subtitle = sprintf("%d samples; without trophic mode (grey): median %.0f%% of the reads per sample",
                            n_distinct(df$run_accession), 100 * median(n_na))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0), panel.spacing.x = unit(1, "pt"),
          legend.position = "right")
}
save_plot(plot_modes(filter(shares, set == "all eukaryotes"), "Trophic modes of the 18S reads (PR2 names matched to Tara W4)"),
          "05_trophic_modes_eukaryotes", width = 16, height = 6)
save_plot(plot_modes(filter(shares, set == "dinoflagellates"), "Trophic modes of the dinoflagellate reads (Dinophyceae, PR2 -> Tara W4)"),
          "05_trophic_modes_dinoflagellates", width = 16, height = 6)
cat("plots written to", plot_dir, "\n")
