# 02: Sample QC
#   1) read depth per sample (which samples fall below min_sample_reads)
#   2) share of prokaryote / plastid / mitochondrion / eukaryote 18S reads per sample
#   3) negative controls: what is in them, and how abundant those ASVs are in
#      the 2024 samples (same sequencing run as the blanks)
#   4) samples without data (not in the count table)

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/trec_biobank/analysis/metab"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))

asv_info    <- read_prepared("asv_info")
counts_long <- read_prepared("counts_long")
samples     <- read_prepared("samples")

env <- samples %>% filter(sample_type == "environmental", in_count_table)


# ---- 1. Read depth -----------------------------------------------------------
p_depth <- env %>%
  mutate(site = coalesce(site, "no metadata")) %>%
  ggplot(aes(x = reorder(run_accession, total_reads), y = total_reads, fill = sequencing_batch)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = min_sample_reads, linetype = "dashed", colour = "red") +
  facet_grid(~ site, scales = "free_x", space = "free_x") +
  scale_y_log10(labels = scales::label_comma()) +
  labs(x = NULL, y = "reads in ASV table (log)", fill = "sequencing batch",
       title = "Read depth per sample", subtitle = paste("dashed: min_sample_reads =", min_sample_reads)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(angle = 90, hjust = 0))
save_plot(p_depth, "02_read_depth", width = 12, height = 5)

env %>% filter(!pass_qc) %>%
  select(run_accession, site, date, size_fraction, total_reads, has_metadata, notes) %>%
  arrange(total_reads) %>% print()


# ---- 2. Read categories per sample -------------------------------------------
cat_long <- env %>%
  select(run_accession, site, date, time_of_day, size_fraction, replicate, starts_with("reads_")) %>%
  pivot_longer(starts_with("reads_"), names_to = "category", names_prefix = "reads_", values_to = "reads") %>%
  group_by(run_accession) %>% mutate(fraction = reads / sum(reads)) %>% ungroup() %>%
  mutate(category = factor(category, levels = names(category_colours)),
         site = coalesce(site, "no metadata"),
         sample = paste(date, time_of_day, replicate, run_accession))

p_cat <- ggplot(cat_long, aes(x = sample, y = fraction, fill = category)) +
  geom_col(width = 0.9) +
  facet_grid(size_fraction ~ site, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = category_colours) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "share of reads", fill = NULL, title = "Read categories per sample") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 0), legend.position = "bottom")
save_plot(p_cat, "02_read_categories", width = 13, height = 8)

cat_summary <- cat_long %>% group_by(site, category) %>%
  summarise(median_share = median(fraction), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = median_share)
write_tsv(cat_summary, file.path(table_dir, "02_category_share_median_per_site.tsv"))
print(cat_summary)


# ---- 3. Negative controls ----------------------------------------------------
controls <- samples %>% filter(sample_type == "control")
controls %>% select(run_accession, control_name, total_reads, n_asvs, in_count_table) %>% print()

# Every ASV in a blank, with its reads there and its mean relative abundance
# and prevalence in the 2024 environmental samples (same sequencing run).
rel <- counts_long %>% group_by(run_accession) %>% mutate(rel = reads / sum(reads)) %>% ungroup()
env24 <- samples %>% filter(pass_qc, sequencing_batch == "2024") %>% pull(run_accession)

control_asvs <- rel %>%
  filter(run_accession %in% controls$run_accession) %>%
  group_by(ASV_ID) %>% summarise(blank_reads = sum(reads), n_blanks = n(), .groups = "drop") %>%
  left_join(rel %>% filter(run_accession %in% env24) %>% group_by(ASV_ID) %>%
              summarise(env24_mean_rel = sum(rel) / length(env24), env24_prevalence = n() / length(env24), .groups = "drop"),
            by = "ASV_ID") %>%
  mutate(across(c(env24_mean_rel, env24_prevalence), ~ replace_na(.x, 0))) %>%
  left_join(asv_info %>% select(ASV_ID, category, lowest_taxon, taxonomy), by = "ASV_ID") %>%
  arrange(desc(blank_reads))
write_tsv(control_asvs, file.path(table_dir, "02_control_asvs.tsv"))
cat("ASVs in blanks:", nrow(control_asvs), " reads:", sum(control_asvs$blank_reads), "\n")
print(head(control_asvs %>% select(-taxonomy), 15))


# ---- 4. Samples without data -------------------------------------------------
missing <- samples %>% filter(!in_count_table) %>%
  select(run_accession, sample_type, site, date, size_fraction, sequencing_batch, flowcell)
write_tsv(missing, file.path(table_dir, "02_samples_not_in_count_table.tsv"), na = "")
cat("samples not in the count table:", nrow(missing), "\n")
print(count(missing, sample_type, site))
