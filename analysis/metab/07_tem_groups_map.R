# 07: Map of Europe with the 18S reads grouped into the major groups of the TEM screen
#
# The TEM annotation uses eight major groups (group_lookup in trec_biobank/analysis/tem/common.R):
# Dinoflagellates, Diatoms, Ciliates, Flagellates, Coccolithophores, Dictyocha, Others,
# Haptophytes. They sit at different taxonomic levels, so every 18S ASV (category
# eukaryote_nuclear) is put into a group by the rules in `group_rules` below: one PR2 rank
# and the names at that rank. The rules are checked from top to bottom and the first match
# wins; ASVs that match no rule are "Others".
# Metazoa (~half of the 18S reads) are not single cells in the TEM screen; with
# exclude_metazoa = TRUE they are left out of the total, so the shares are shares of the
# non-metazoan 18S reads.
#
# One stacked bar per site x size fraction = mean share over the samples of that site and
# size fraction (every sample weighted equally). Map layout: map_common.R.
# Input:  data/ of 01_prepare_data.R
# Output: plots/07_tem_groups_map
#         tables/07_tem_group_per_sample.tsv     share per sample and group
#         tables/07_tem_groups_site_fraction.tsv mean share per site, size fraction and group
#         tables/07_tem_group_lineages.tsv       which PR2 lineages went into which group (to check the rules)

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/trec_biobank/analysis/metab"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "map_common.R"))

# ---- settings ----------------------------------------------------------------
exclude_metazoa <- TRUE

# PR2 5.1 names. rank = column pr2_<rank> of asv_info. First match wins, so the
# narrow groups (Dictyocha, Coccolithophores) come before their parents.
group_rules <- tribble(
  ~group,             ~rank,          ~names,
  "Metazoa",          "Subdivision",  "Metazoa",
  "Dictyocha",        "Genus",        c("Dictyocha", "Dictyochales_X"),
  # calcifying haptophytes: the coccolithophore orders, plus Emiliania / Gephyrocapsa
  # (Noelaerhabdaceae, Isochrysidales) and Braarudosphaera
  "Coccolithophores", "Order",        c("Coccolithales", "Syracosphaerales", "Zygodiscales", "Calcihaptophycidae"),
  "Coccolithophores", "Family",       c("Noelaerhabdaceae", "Braarudosphaeraceae"),
  "Haptophytes",      "Division",     "Haptophyta",
  # incl. Syndiniales (parasitic) and Noctilucophyceae
  "Dinoflagellates",  "Subdivision",  "Dinoflagellata",
  "Diatoms",          "Class",        c("Mediophyceae", "Coscinodiscophyceae", "Bacillariophyceae", "Diatomeae_X"),
  "Ciliates",         "Subdivision",  "Ciliophora",
  # flagellates: green algae, cryptophytes, the flagellate stramenopile classes and
  # heterotrophic flagellate lineages (MAST / bicosoecids, choanoflagellates, euglenozoa, ...)
  "Flagellates",      "Division",     c("Chlorophyta", "Prasinodermophyta", "Cryptophyta", "Kathablepharidacea",
                                        "Telonemia", "Picozoa", "Apusomonada"),
  "Flagellates",      "Subdivision",  c("Bigyra", "Choanoflagellata", "Euglenozoa"),
  "Flagellates",      "Class",        c("Chrysophyceae", "Pelagophyceae", "Dictyochophyceae", "Bolidophyceae",
                                        "Pinguiophyceae", "Raphidophyceae", "MOCH-1", "MOCH-2", "MOCH-3",
                                        "MOCH-4", "MOCH-5")
)
# order of the TEM screen, Haptophytes next to Coccolithophores and Others last
group_levels  <- c("Dinoflagellates", "Diatoms", "Ciliates", "Flagellates", "Coccolithophores",
                   "Haptophytes", "Dictyocha", "Others")
group_colours <- c(Dinoflagellates = "#C44E52", Diatoms = "#8C6D31", Ciliates = "#4C72B0",
                   Flagellates = "#55A868", Coccolithophores = "#8172B2", Haptophytes = "#CCB9E8",
                   Dictyocha = "#E1A100", Others = "grey75")

# ---- data --------------------------------------------------------------------
asv_info    <- read_prepared("asv_info") %>% filter(category == "eukaryote_nuclear")
counts_long <- read_prepared("counts_long")
samples     <- read_prepared("samples") %>% filter(pass_qc)

## group per ASV: apply the rules from top to bottom
asv_group <- asv_info %>% mutate(group = NA_character_)
for (i in seq_len(nrow(group_rules))) {
  col <- paste0("pr2_", group_rules$rank[i])
  hit <- is.na(asv_group$group) & asv_group[[col]] %in% group_rules$names[[i]]
  asv_group$group[hit] <- group_rules$group[i]
}
asv_group <- asv_group %>% mutate(group = coalesce(group, "Others"))
if (exclude_metazoa) asv_group <- asv_group %>% filter(group != "Metazoa")
stopifnot(all(asv_group$group %in% group_levels))

counts <- counts_long %>%
  filter(run_accession %in% samples$run_accession) %>%
  inner_join(asv_group %>% select(ASV_ID, group), by = "ASV_ID")

## lineages per group, to check the rules
lineages <- counts %>%
  count(ASV_ID, wt = reads, name = "reads") %>%
  left_join(asv_group, by = "ASV_ID") %>%
  group_by(group, pr2_Division, pr2_Subdivision, pr2_Class, pr2_Order) %>%
  summarise(n_asvs = n(), reads = sum(reads), .groups = "drop") %>%
  mutate(share = reads / sum(reads)) %>%
  arrange(factor(group, group_levels), desc(reads))
write_tsv(lineages, file.path(table_dir, "07_tem_group_lineages.tsv"), na = "")

cat("\nshare of the", if (exclude_metazoa) "non-metazoan", "18S reads per group (all samples):\n")
print(lineages %>% count(group, wt = reads, name = "reads") %>%
        mutate(share = round(reads / sum(reads), 3)) %>% arrange(factor(group, group_levels)))

# ---- share per sample and per site x size fraction ---------------------------
per_sample <- counts %>%
  count(run_accession, group, wt = reads, name = "reads") %>%
  complete(run_accession, group = group_levels, fill = list(reads = 0)) %>%
  group_by(run_accession) %>% mutate(total_reads = sum(reads), share = reads / total_reads) %>% ungroup() %>%
  filter(total_reads > 0) %>%
  left_join(samples %>% select(run_accession, site, date, time_of_day, size_fraction, replicate, label),
            by = "run_accession")
write_tsv(per_sample, file.path(table_dir, "07_tem_group_per_sample.tsv"), na = "")

site_fraction <- per_sample %>%
  filter(!is.na(site), !is.na(size_fraction)) %>%
  group_by(site, size_fraction, group) %>%
  summarise(share = mean(share), n_samples = n_distinct(run_accession), .groups = "drop") %>%
  mutate(group = factor(group, group_levels))
write_tsv(site_fraction, file.path(table_dir, "07_tem_groups_site_fraction.tsv"))

# ---- map ---------------------------------------------------------------------
check_site_coords(site_fraction$site)
p <- bar_map(site_fraction, fill = "group", colours = group_colours, legend_title = "TEM major group",
             title = paste("18S reads in the major groups of the TEM screen",
                           if (exclude_metazoa) "(without Metazoa)"),
             subtitle = paste("one bar per size fraction (um; number of samples): mean share of the",
                              if (exclude_metazoa) "non-metazoan", "18S reads over the samples"),
             legend_nrow = 1)
save_map(p, "07_tem_groups_map")
cat("plots written to", plot_dir, "\n")
