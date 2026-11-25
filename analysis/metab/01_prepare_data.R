# 01: Read the taxonomy workflow output and the sample mapping once, and save
# three tidy objects that all later scripts use:
#
#   data/asv_info.rds     one row per ASV: taxonomy, category, plotting groups
#   data/counts_long.rds  one row per ASV x sample with reads > 0
#   data/samples.rds      one row per run: metadata + read totals per category
#
# Run this again whenever taxonomy_dir or sample_mapping.tsv changes.

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/trec_biobank/analysis/metab"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))

# ---- 1. ASV taxonomy ---------------------------------------------------------
asv_tax <- read_tsv(file.path(taxonomy_dir, "asv_taxonomy.tsv"),
                    col_types = cols(.default = col_character()), na = c("", "NA"))

# Columns of the organelle step exist only in runs made after it was added.
has_host <- "organelle_host_taxonomy" %in% names(asv_tax)
cat("organelle host columns present:", has_host, "\n")
if (!has_host) asv_tax <- asv_tax %>% mutate(organelle_host_taxonomy = NA_character_, organelle_host_lowest_taxon = NA_character_)

# Host lineage split into PR2 ranks (":plas"/":mito" tags are already removed).
pr2_ranks <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
host <- str_split_fixed(coalesce(asv_tax$organelle_host_taxonomy, ""), ";", length(pr2_ranks))
host[host == ""] <- NA
colnames(host) <- paste0("host_", pr2_ranks)

asv_info <- asv_tax %>%
  select(ASV_ID, ssu, category, reference, taxonomy, lowest_rank, lowest_taxon, length,
         starts_with("silva_") & !starts_with("silva_boot"), starts_with("pr2_") & !starts_with("pr2_boot"),
         organelle_host_taxonomy, organelle_host_lowest_taxon) %>%
  bind_cols(as_tibble(host)) %>%
  mutate(
    length = as.integer(length),
    # One "plotting group" per ASV at a coarse and a finer level, from the
    # reference that fits its category.
    group_coarse = case_when(
      # PR2 Subdivision: Metazoa, Fungi, Dinoflagellata, Ciliophora, Gyrista (diatoms), ...
      category == "eukaryote_nuclear" ~ coalesce(pr2_Subdivision, pr2_Division, pr2_Supergroup, "unclassified eukaryote"),
      category == "prokaryote"        ~ coalesce(silva_Phylum, silva_Kingdom, "unclassified prokaryote"),
      category == "plastid"           ~ coalesce(host_Subdivision, host_Division, "plastid, host unknown"),
      category == "mitochondrion"     ~ coalesce(host_Subdivision, host_Division, "mitochondrion, host unknown"),
      TRUE                            ~ "unassigned"),
    group_fine = case_when(
      category == "eukaryote_nuclear" ~ coalesce(pr2_Class, paste(group_coarse, "(class unknown)")),
      category == "prokaryote"        ~ coalesce(silva_Class, paste(group_coarse, "(class unknown)")),
      category %in% c("plastid", "mitochondrion") ~ coalesce(host_Class, paste(group_coarse, "(class unknown)")),
      TRUE                            ~ "unassigned")
  )

asv_info %>% count(category)


# ---- 2. Count table -> long format -------------------------------------------
mapping <- read_tsv(mapping_f, col_types = cols(.default = col_character()))

counts_wide <- data.table::fread(file.path(taxonomy_dir, "asv_table_taxonomy.tsv"),
                                 select = c("ASV_ID", intersect(
                                   names(data.table::fread(file.path(taxonomy_dir, "asv_table_taxonomy.tsv"), nrows = 0)),
                                   mapping$run_accession)))
cat("ASVs:", nrow(counts_wide), " samples in count table:", ncol(counts_wide) - 1, "\n")

counts_long <- data.table::melt(counts_wide, id.vars = "ASV_ID", variable.name = "run_accession",
                                value.name = "reads", variable.factor = FALSE)
counts_long <- as_tibble(counts_long[reads > 0])
rm(counts_wide); invisible(gc())
stopifnot(all(counts_long$ASV_ID %in% asv_info$ASV_ID))


# ---- 3. Samples: metadata + read totals --------------------------------------
per_category <- counts_long %>%
  left_join(asv_info %>% select(ASV_ID, category), by = "ASV_ID") %>%
  count(run_accession, category, wt = reads) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0, names_prefix = "reads_")

per_sample <- counts_long %>%
  group_by(run_accession) %>%
  summarise(total_reads = sum(reads), n_asvs = n(), .groups = "drop")

samples <- mapping %>%
  mutate(date = as.Date(date), tara = as.logical(tara)) %>%
  select(-any_of(c("in_asv_table", "asv_reads"))) %>%        # recomputed here from the current run
  left_join(per_sample, by = "run_accession") %>%
  left_join(per_category, by = "run_accession") %>%
  mutate(across(c(total_reads, n_asvs, starts_with("reads_")), ~ replace_na(.x, 0)),
         in_count_table = run_accession %in% counts_long$run_accession,
         year = format(date, "%Y"),
         # A readable label for plots, e.g. "Bilbao 2023-10-01 AM 10-40 R1"
         label = if_else(sample_type == "control", paste("control", control_name, run_accession),
                         paste(site, date, time_of_day, size_fraction, replicate)),
         pass_qc = sample_type == "environmental" & has_metadata == "TRUE" & total_reads >= min_sample_reads)

samples %>% count(sample_type, in_count_table, pass_qc)


# ---- 4. Save -----------------------------------------------------------------
saveRDS(asv_info,    file.path(data_dir, "asv_info.rds"))
saveRDS(counts_long, file.path(data_dir, "counts_long.rds"))
saveRDS(samples,     file.path(data_dir, "samples.rds"))
write_tsv(samples, file.path(table_dir, "samples_with_read_totals.tsv"), na = "")
cat("saved to", data_dir, "\n")
