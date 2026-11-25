## Map of Europe: imaging-based (TEM) and sequencing-based (metabarcoding) trophic groups per site
## per site a group of stacked bars on the site's position:
##   TEM    (black outline)  chloroplast + food vacuole / chloroplast only / no chloroplast
##                           of the dinoflagellate cells (11_chloro_foodvac_sites.R)
##   metaB  one bar per size fraction (metab_fractions): mean share of the trophic modes
##          (Tara W4) of the dinoflagellate reads over the samples
##          (analysis/metab/05_trophic_modes.R + 06_trophic_modes_map.R)
## sites without TEM cells only get the metaB bars. Sites too close to each other are moved
## by dx / dy (m) with a line to the real position.
## Input:  tables/chloro_foodvac_site_counts.tsv of 11 (this run)
##         metab table 06_trophic_modes_site_fraction.tsv (env TEM_METAB_TROPHIC_TABLE)
## Output: plots/trophic_map_combined               TEM: one bar per site (all size fractions)
##         plots/trophic_map_combined_by_fraction   TEM and metaB both per size fraction (pairs of bars
##             TEM | metaB per size fraction), site labels with the sampling month, and on the right
##             the sampling dates of every site (one panel per site, own day axis) (TEM: dates of the annotated cells, metaB: dates
##             of the samples of the shown size fractions; metab table 05_trophic_mode_per_sample.tsv,
##             env TEM_METAB_SAMPLE_TABLE)
##         tables/trophic_map_combined.tsv, trophic_map_combined_by_fraction.tsv, trophic_map_sampling_dates.tsv

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- NA_character_   # not read here (only the tables of 11 and of the metab analysis)
default_outdir     <- "C:/projects/trec_biobank/analysis/tem/rstudio"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/trec_biobank/analysis/tem"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(patchwork)
})
source(file.path(script_dir, "map_common.R"))   # site_coords (bar shifts dx / dy), map_with_bars(), timeline_plot()
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)
## default: the newest metab run (<data_root>/analysis/metab/<run>/) that has the table (paths.R)
default_metab_table <- file.path(trec_latest_run(trec_paths$metab_runs, "tables/06_trophic_modes_site_fraction.tsv"),
                                 "tables/06_trophic_modes_site_fraction.tsv")
metab_table <- trec_env("TEM_METAB_TROPHIC_TABLE", default_metab_table)
metab_sample_table <- Sys.getenv("TEM_METAB_SAMPLE_TABLE", file.path(dirname(metab_table), "05_trophic_mode_per_sample.tsv"))

## ---- settings ----
## (site coordinates, bar shifts and bar sizes: map_common.R)
metab_fractions <- c("5-40", "10-40", "10-200", "20-200")   # metaB size fractions shown
metab_set       <- "dinoflagellates"                        # or "all eukaryotes"

tem_levels <- c("chloroplast + food vacuole", "chloroplast only", "no chloroplast")
## aligned with metaB: both chloroplast groups green (= phototroph, light / dark), no chloroplast = heterotroph orange
tem_cols   <- c(`chloroplast + food vacuole` = "#7CC8B0", `chloroplast only` = "#1B9E77", `no chloroplast` = "#D95F02")
metab_levels <- c("phototroph", "photosymbiosis", "heterotroph", "Metazoa", "NA: genus not in W4", "NA: no genus (PR2)")
metab_cols   <- c(phototroph = "#1B9E77", photosymbiosis = "#A6D854", heterotroph = "#D95F02",
                  Metazoa = "#5B7DB1", `NA: genus not in W4` = "grey55", `NA: no genus (PR2)` = "grey82")

## ---- data ----
tem_file  <- file.path(opt$table_dir, "chloro_foodvac_site_counts.tsv")
cell_file <- file.path(opt$table_dir, "chloro_foodvac_cells.tsv")
for (f in c(tem_file, cell_file)) if (!file.exists(f)) stop("Missing ", f, ": run 11_chloro_foodvac_sites.R first")
for (f in c(metab_table, metab_sample_table)) if (!file.exists(f))
  stop("Missing ", f, ": run analysis/metab steps 05 + 06 first (or set TEM_METAB_TROPHIC_TABLE / TEM_METAB_SAMPLE_TABLE)")

tem_site <- read_tsv(tem_file, show_col_types = FALSE) %>%
  transmute(site = tem_site_name(site), source = "TEM", share = frac, frac = NA_character_, n, cat = troph) %>%
  group_by(site) %>% mutate(n = sum(n)) %>% ungroup()
tem_cells <- read_tsv(cell_file, show_col_types = FALSE) %>%
  mutate(site = tem_site_name(site), frac = frac_label(size_frac))
tem_frac <- tem_cells %>%
  count(site, frac, troph) %>%
  group_by(site, frac) %>% mutate(n_cells = sum(n), share = n / n_cells) %>% ungroup() %>%
  transmute(site, source = "TEM", frac, n = n_cells, cat = troph, share)

metab <- read_tsv(metab_table, show_col_types = FALSE) %>%
  filter(set == metab_set, size_fraction %in% metab_fractions)
unknown <- setdiff(unique(metab$site), site_coords$site)
if (length(unknown) > 0) stop("metaB site(s) without coordinates: ", paste(unknown, collapse = ", "))
metab <- metab %>% transmute(site, source = "metaB", frac = size_fraction, n = n_samples, cat = mode, share)

## sampling dates: TEM = dates of the annotated cells, metaB = dates of the samples shown
dates <- bind_rows(
  tem_cells %>% distinct(site, date = as.Date(sampling_date)) %>% mutate(source = "TEM"),
  read_tsv(metab_sample_table, show_col_types = FALSE) %>%
    filter(set == metab_set, sample_type == "environmental", !is.na(site), size_fraction %in% metab_fractions) %>%
    distinct(site, date = as.Date(date)) %>% mutate(source = "metaB")) %>%
  filter(!is.na(date))
write_tsv(dates %>% arrange(site, source, date), file.path(opt$table_dir, "trophic_map_sampling_dates.tsv"))

## two legends (TEM groups / metaB trophic modes)
legends_of <- function(bars) {
  wrap_elements(legend_of(bars, tem_levels, tem_cols, "TEM: dinoflagellate cells", nrow = 3)) |
    wrap_elements(legend_of(bars, metab_levels, metab_cols, "metaB: trophic mode of the reads (Tara W4)", nrow = 2))
}
metab_reads <- if (metab_set == "dinoflagellates") "dinoflagellate" else "18S"

## ---- 1) TEM per site (all size fractions) + metaB per size fraction ----
bars1 <- bind_rows(tem_site, metab) %>%
  mutate(cat = factor(cat, c(tem_levels, metab_levels)),
         bar = coalesce(frac, "TEM"), pair = source,
         bar_order = if_else(source == "TEM", 0L, match(frac, metab_fractions)),
         label = if_else(source == "TEM", sprintf("TEM (%d cells)", n), sprintf("%s (%d)", frac, n)))
m1 <- map_with_bars(
  bars1, "Trophic groups of dinoflagellates: imaging (TEM) and metabarcoding (metaB)",
  sprintf(paste0("per site: TEM bar (black outline; share of the annotated cells), then one metaB bar per size fraction\n",
                 "(%s um; number of samples; mean share of the %s reads, trophic modes from Tara W4)"),
          paste(metab_fractions, collapse = ", "), metab_reads),
  site_labels = setNames(site_coords$site, site_coords$site), fill_values = c(tem_cols, metab_cols))
p1 <- m1$plot / legends_of(m1$bars) + plot_layout(heights = c(1, 0.08)) & white_bg
save_plot(p1, "trophic_map_combined", 11, 11 * m1$asp + 1.6, opt$plot_dir)
write_tsv(bars1 %>% select(site, source, size_fraction = frac, n, cat, share), file.path(opt$table_dir, "trophic_map_combined.tsv"))

## ---- 2) TEM and metaB per size fraction + sampling time span ----
site_months <- dates %>% group_by(site) %>% summarise(m = fmt_months(date), .groups = "drop")
site_labels2 <- setNames(paste0(site_months$site, "\n", site_months$m), site_months$site)

bars2 <- bind_rows(tem_frac, metab) %>%
  mutate(cat = factor(cat, c(tem_levels, metab_levels)),
         bar = paste(source, frac), pair = frac,
         bar_order = frac_order(frac) * 2L + if_else(source == "TEM", 0L, 1L),
         label = if_else(source == "TEM", sprintf("TEM %s (%d)", frac, n), sprintf("metaB %s (%d)", frac, n)))
m2 <- map_with_bars(
  bars2, "Trophic groups of dinoflagellates per size fraction: imaging (TEM) and metabarcoding (metaB)",
  sprintf(paste0("per site and size fraction a pair of bars: TEM (black outline; number of cells) | metaB (number of samples; ",
                 "mean share of the %s reads, Tara W4)\nmetaB size fractions: %s um; TEM: all annotated size fractions"),
          metab_reads, paste(metab_fractions, collapse = ", ")),
  site_labels = site_labels2, fill_values = c(tem_cols, metab_cols), w = bar_w_frac, gap = bar_gap_frac, bar_h = bar_h_frac)
p_time <- timeline_plot(dates)

p2 <- (m2$plot | p_time) + plot_layout(widths = c(3, 1))
p2 <- p2 / legends_of(m2$bars) + plot_layout(heights = c(1, 0.08)) & white_bg
## map panel width ~ (figure width - axis labels of the timeline) * 3/4
save_plot(p2, "trophic_map_combined_by_fraction", 16, (16 - 1) * 3 / 4 * m2$asp + 2.4, opt$plot_dir)
write_tsv(bars2 %>% select(site, source, size_fraction = frac, n, cat, share),
          file.path(opt$table_dir, "trophic_map_combined_by_fraction.tsv"))
