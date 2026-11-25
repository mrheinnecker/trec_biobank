## TEM only: maps of Europe with organelle groups of the dinoflagellate cells per site and size fraction
## (same layout as trophic_map_combined_by_fraction of 13, without metabarcoding)
##   1) crystalline reticulate network: cryst_rich_RN present / absent
##   2) storage: starch + lipid droplets / starch only / lipid droplets only / none
## one stacked bar per site x size fraction (share of the cells; number of cells under the bar),
## site labels with the sampling month, sampling dates per site on the right.
## Cells with a missing value in a grouping column are left out of that map (uncertain "?" = 0, as everywhere).
## Output: plots/tem_map_cryst_rich_RN, plots/tem_map_starch_lipid
##         tables/tem_map_groups_cells.tsv (group per cell), tables/tem_map_groups_counts.tsv

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- NA_character_   # NA = read the main sheet from Google Sheets; or a tsv snapshot
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
## the bar shifts in map_common.R make room for the metaB bars; with TEM only all bars stand on their site
sites_all <- sites_all %>% mutate(dx = 0, dy = 0, bx = X, by = Y)

## ---- settings ----
groupings <- list(
  cryst_rich_RN = list(
    title  = "Crystalline reticulate network (cryst_rich_RN) in dinoflagellate cells (TEM)",
    legend = "TEM: crystalline reticulate network",
    levels = c("crystalline RN", "no crystalline RN"),
    cols   = c(`crystalline RN` = "#6A51A3", `no crystalline RN` = "#CBC9E2"),
    group  = function(d) if_else(d$cryst_rich_RN == 1, "crystalline RN", "no crystalline RN"),
    needs  = "cryst_rich_RN"),
  starch_lipid = list(
    title  = "Starch and lipid droplets in dinoflagellate cells (TEM)",
    legend = "TEM: storage products",
    levels = c("starch + lipid droplets", "starch only", "lipid droplets only", "none"),
    cols   = c(`starch + lipid droplets` = "#8C510A", `starch only` = "#DFC27D",
               `lipid droplets only` = "#5AB4AC", none = "#E3E3E3"),
    group  = function(d) case_when(d$starch == 1 & d$lipid_drop == 1 ~ "starch + lipid droplets",
                                   d$starch == 1                      ~ "starch only",
                                   d$lipid_drop == 1                  ~ "lipid droplets only",
                                   TRUE                               ~ "none"),
    needs  = c("starch", "lipid_drop"))
)

## ---- data ----
dino <- load_dino_data(opt) %>%
  mutate(site_name = tem_site_name(site), frac = frac_label(size_frac),
         date = as.Date(str_match(source_name, "_\\d{9}_[A-Z]{3}_[^_]+_(\\d{8})_")[, 2], "%Y%m%d"))

cells <- imap_dfr(groupings, function(g, key) {
  d <- dino %>% filter(if_all(all_of(g$needs), ~ !is.na(.x)))
  tibble(grouping = key, name = d$name, site = d$site_name, frac = d$frac, date = d$date, group = g$group(d))
})
counts <- cells %>%
  count(grouping, site, frac, group) %>%
  group_by(grouping, site, frac) %>% mutate(n_cells = sum(n), share = n / n_cells) %>% ungroup()
write_tsv(cells, file.path(opt$table_dir, "tem_map_groups_cells.tsv"))
write_tsv(counts, file.path(opt$table_dir, "tem_map_groups_counts.tsv"))

## site labels with the sampling month(s), sampling dates of the annotated cells
dates <- dino %>% filter(!is.na(date)) %>% distinct(site = site_name, date) %>% mutate(source = "TEM")
site_months <- dates %>% group_by(site) %>% summarise(m = fmt_months(date), .groups = "drop")
site_labels <- setNames(paste0(site_months$site, "\n", site_months$m), site_months$site)
p_time <- timeline_plot(dates, subtitle_sources = "dates of the annotated cells")

## ---- one map per grouping ----
for (key in names(groupings)) {
  g <- groupings[[key]]
  bars <- counts %>%
    filter(grouping == key) %>%
    transmute(site, source = "cells", frac, n = n_cells, cat = factor(group, g$levels), share,
              bar = frac, pair = frac, bar_order = frac_order(frac), label = sprintf("%s (%d)", frac, n_cells))
  m <- map_with_bars(
    bars, g$title,
    sprintf("one bar per site and size fraction (um; number of cells): share of the annotated dinoflagellate cells; %d cells",
            n_distinct(cells$name[cells$grouping == key])),
    site_labels = site_labels, fill_values = g$cols, w = bar_w_frac, gap = bar_gap_frac, bar_h = bar_h_frac)
  legend <- wrap_elements(legend_of(m$bars, g$levels, g$cols, g$legend, nrow = 1))
  p <- (m$plot | p_time) + plot_layout(widths = c(3, 1))
  p <- p / legend + plot_layout(heights = c(1, 0.05)) & white_bg
  save_plot(p, paste0("tem_map_", key), 16, (16 - 1) * 3 / 4 * m$asp + 2, opt$plot_dir)
}
