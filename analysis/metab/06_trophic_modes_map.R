# 06: Map of Europe with the trophic modes of 05_trophic_modes.R per site and size fraction
#
# One stacked bar per site x size fraction = mean share of the trophic modes over the samples
# of that site and size fraction (every sample weighted equally). Map, site coordinates and
# bar layout: map_common.R.
# Input:  tables/05_trophic_mode_per_sample.tsv of 05_trophic_modes.R (run 05 first)
# Output: plots/06_trophic_modes_map_dinoflagellates, plots/06_trophic_modes_map_eukaryotes,
#         tables/06_trophic_modes_site_fraction.tsv

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/trec_biobank/analysis/metab"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "map_common.R"))

# ---- settings ----------------------------------------------------------------
# from 05_trophic_modes.R
mode_levels <- c("phototroph", "photosymbiosis", "heterotroph", "Metazoa",
                 "NA: genus not in W4", "NA: no genus (PR2)")
mode_colours <- c(phototroph = "#1B9E77", photosymbiosis = "#A6D854", heterotroph = "#D95F02",
                  Metazoa = "#5B7DB1", `NA: genus not in W4` = "grey55", `NA: no genus (PR2)` = "grey82")

# ---- data --------------------------------------------------------------------
shares_f <- file.path(table_dir, "05_trophic_mode_per_sample.tsv")
if (!file.exists(shares_f)) stop("Missing ", shares_f, ": run 05_trophic_modes.R first")
shares <- read_tsv(shares_f, show_col_types = FALSE) %>%
  filter(sample_type == "environmental", !is.na(site), !is.na(size_fraction))

check_site_coords(shares$site)

## mean share per site x size fraction x mode; size fractions ordered by lower, then upper bound
site_fraction <- shares %>%
  group_by(set, site, size_fraction, mode) %>%
  summarise(share = mean(share), n_samples = n_distinct(run_accession), .groups = "drop") %>%
  mutate(mode = factor(mode, mode_levels))
write_tsv(site_fraction, file.path(table_dir, "06_trophic_modes_site_fraction.tsv"))

# ---- plot --------------------------------------------------------------------
world <- base_map()
map_plot <- function(set_name, title) {
  bar_map(filter(site_fraction, set == set_name), fill = "mode", colours = mode_colours,
          legend_title = "trophic mode (W4)", title = title, world = world)
}

save_map(map_plot("dinoflagellates", "Trophic modes of the dinoflagellate reads (PR2 -> Tara W4)"),
         "06_trophic_modes_map_dinoflagellates")
save_map(map_plot("all eukaryotes", "Trophic modes of the 18S reads (PR2 -> Tara W4)"),
         "06_trophic_modes_map_eukaryotes")
cat("plots written to", plot_dir, "\n")
