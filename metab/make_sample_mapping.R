#!/usr/bin/env Rscript
# Main sample mapping table for the TREC metabarcoding data (16S V4V5).
#
# Joins the two raw tables:
#   raw_input/sampleinfo.tsv     ENA file list: one row per FASTQ file (R1 + R2),
#                                with run accession (ERR...), BioSample (SAMEA...),
#                                dataset, sample type and submitted filename.
#   raw_input/metab_samples.tsv  sampling sheet: one row per sample, with
#                                Barcode ID, site, date, time, size fraction, ...
#
# The key: the Barcode ID in the sampling sheet is the BioSample number without
# the "SAMEA" prefix (Barcode 112518358 = SAMEA112518358 = TREC112518358).
#
# Output: sample_mapping.tsv, one row per sequencing run (= one ampliseq sample).
#
# Written to be run step by step in RStudio; it also runs with Rscript.
# Paths are set in section 0.

library(tidyverse)

# ---- 0. Paths ----------------------------------------------------------------
base_dir      <- "C:/projects/trec_biobank/metab"
sampleinfo_f  <- file.path(base_dir, "raw_input/sampleinfo.tsv")
samples_f     <- file.path(base_dir, "raw_input/metab_samples.tsv")
# optional: ampliseq/taxonomy count table, used to add the reads per run
asv_table_f   <- file.path(base_dir, "ampliseq_2026-11-03/results/taxonomy/asv_table_taxonomy.tsv")
out_f         <- file.path(base_dir, "sample_mapping.tsv")


# ---- 1. ENA file list: one row per FASTQ file -> one row per run -------------
# Everything is read as character so IDs and dates stay unchanged.
ena_files <- read_tsv(sampleinfo_f, col_types = cols(.default = col_character()))

# The Description field holds "key value; key value; ...". Pull out the parts we
# need. The 2023 files have no "library construction protocol" entry.
ena_files <- ena_files %>%
  mutate(
    run_accession      = str_extract(Name, "^ERR[0-9]+"),
    read               = str_extract(Name, "(?<=_)[12](?=\\.fastq)"),
    study_accession    = str_extract(Description, "(?<=study accession )[^;]+"),
    library_protocol   = str_extract(Description, "(?<=library construction protocol )[^;]+"),
    submitted_filename = str_extract(Description, "(?<=submitted filename )[^;]+")
  )

# Check: every run has exactly one R1 and one R2 file.
ena_files %>% count(run_accession, read) %>% count(read)
stopifnot(all(count(ena_files, run_accession)$n == 2))

# Reshape to one row per run: put R1/R2 file name and MD5 side by side.
ena_runs <- ena_files %>%
  select(run_accession, read, Name, MD5, biosample = Samples, sample_type = type,
         ena_dataset = Analyses, study_accession, library_protocol, submitted_filename) %>%
  pivot_wider(names_from = read, values_from = c(Name, MD5, submitted_filename),
              names_glue = "{.value}_R{read}") %>%
  rename(fastq_R1 = Name_R1, fastq_R2 = Name_R2, md5_R1 = MD5_R1, md5_R2 = MD5_R2)

# The submitted filename tells us more, e.g.
#   SAMEA112518358_METAB-16S-V4V5_HGH7NDRX5.UDI229-BID05-2_R1_clean.fastq.gz
#   Tneg_Tag_eau_Ambion_METAB-16S-V4V5_HCJF5DRX7.UDI140-BID18-1_R1.fastq.gz
#   ^ sample / control name      ^ flowcell ^ index/barcode tag
ena_runs <- ena_runs %>%
  mutate(
    submitted_name = str_extract(submitted_filename_R1, "^.+?(?=_METAB-)"),
    flowcell       = str_extract(submitted_filename_R1, "(?<=METAB-16S-V4V5_)[A-Z0-9]+(?=\\.)"),
    index_tag      = str_extract(submitted_filename_R1, "(?<=\\.)UDI[0-9]+-BID[0-9]+-[0-9]+"),
    # 2023: the first ENA submission; 2024: "..._2024" and "analysis-1"
    sequencing_batch = if_else(ena_dataset == "MetaB_16SV4V5_Plankton_AML40um", "2023", "2024"),
    control_name   = if_else(sample_type == "control", submitted_name, NA_character_),
    barcode_id     = str_remove(biosample, "^SAMEA")
  ) %>%
  select(-submitted_filename_R2)

# Checks: one row per run, one run per BioSample, and for real samples the
# submitted filename starts with the BioSample.
stopifnot(!anyDuplicated(ena_runs$run_accession), !anyDuplicated(ena_runs$biosample))
ena_runs %>% filter(sample_type != "control", submitted_name != biosample)   # expect 0 rows
ena_runs %>% count(sequencing_batch, flowcell, sample_type)


# ---- 2. Sampling sheet: one row per sample -----------------------------------
samples_raw <- read_tsv(samples_f, col_types = cols(.default = col_character()))

# The sheet has empty spacer rows, a note with a link at the bottom, and one
# COPAS-sorted Protoperidinium sample without a Barcode ID. Only rows with a
# Barcode ID can be linked to ENA.
samples_raw %>% filter(is.na(`Barcode ID`)) %>% select(`Bottle ID/Event`, `LabID name`, Site, Date, Fraction)

samples <- samples_raw %>%
  filter(!is.na(`Barcode ID`)) %>%
  transmute(
    barcode_id        = `Barcode ID`,
    lab_id            = `LabID name`,
    site              = Site,
    date              = as.Date(Date, format = "%Y%m%d"),
    time_of_day       = Time,                    # AM / MID / PM
    size_fraction_raw = Fraction,
    # "10to40" and "10-40" mean the same; write all as "10-40" (µm)
    size_fraction     = str_replace(Fraction, "to", "-"),
    tara              = as.logical(TARA),
    replicate         = Replicate,
    depth_m           = `Depth (m)`,             # free text (e.g. "5-7", "0-30 m")
    sent_to_genoscope = as.logical(`Sent to Genoscope`),
    notes             = Notes
  )

stopifnot(!anyDuplicated(samples$barcode_id))
stopifnot(!any(is.na(samples$date)))
samples %>% count(size_fraction_raw, size_fraction)


# ---- 3. Join runs and samples ------------------------------------------------
mapping <- ena_runs %>%
  left_join(samples, by = "barcode_id") %>%
  mutate(trec_name    = if_else(sample_type == "control", NA_character_, paste0("TREC", barcode_id)),
         has_metadata = !is.na(site))

# What did not match?
# a) runs without sampling metadata: the 6 negative controls (expected) and
#    2 environmental runs whose barcode is missing from the sampling sheet.
mapping %>% filter(!has_metadata) %>% select(run_accession, biosample, sample_type, ena_dataset, submitted_name)
# b) samples in the sheet without a run on ENA (mostly "low reads" /
#    "to be re-sequenced" / not received by Genoscope).
samples %>% anti_join(ena_runs, by = "barcode_id") %>% count(site, notes)


# ---- 4. Optional: total reads per run in the final ASV table -----------------
# Run IDs are the sample columns of asv_table_taxonomy.tsv. Runs that ampliseq
# dropped (too few reads after filtering) get 0.
if (file.exists(asv_table_f)) {
  asv_header <- names(read_tsv(asv_table_f, n_max = 0, show_col_types = FALSE))
  run_cols   <- intersect(asv_header, mapping$run_accession)
  asv_counts <- read_tsv(asv_table_f, col_select = all_of(run_cols),
                         col_types = cols(.default = col_double()))
  reads <- tibble(run_accession = run_cols, asv_reads = colSums(asv_counts))
  mapping <- mapping %>%
    left_join(reads, by = "run_accession") %>%
    mutate(in_asv_table = run_accession %in% run_cols,
           asv_reads    = replace_na(asv_reads, 0))
  mapping %>% count(sample_type, in_asv_table)
}


# ---- 5. Order columns and write ----------------------------------------------
mapping <- mapping %>%
  select(run_accession, biosample, barcode_id, trec_name, lab_id, sample_type, control_name,
         site, date, time_of_day, size_fraction, replicate, depth_m, tara, notes,
         sequencing_batch, flowcell, index_tag, ena_dataset, study_accession, library_protocol,
         any_of(c("in_asv_table", "asv_reads")), has_metadata,
         size_fraction_raw, sent_to_genoscope, submitted_name, submitted_filename_R1,
         fastq_R1, fastq_R2, md5_R1, md5_R2) %>%
  arrange(sample_type, site, date, time_of_day, size_fraction, replicate)

write_tsv(mapping, out_f, na = "")
cat("wrote", nrow(mapping), "runs to", out_f, "\n")

# Quick overview: runs per site and year
mapping %>% count(sample_type, site, year = format(date, "%Y")) %>% print(n = Inf)
