## Map of Europe: imaging-based (TEM) and sequencing-based (metabarcoding) trophic groups per site
## per site a group of stacked bars on the site's position:
##   TEM    (black outline)  chloroplast + food vacuole / chloroplast only / no chloroplast
##                           of the dinoflagellate cells (11_chloro_foodvac_sites.R)
##   metaB  one bar per size fraction (metab_fractions): mean share of the trophic modes
##          (Tara W4) of the dinoflagellate reads over the samples
##          (metab/analysis/05_trophic_modes.R + 06_trophic_modes_map.R)
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
default_outdir     <- "C:/projects/tem_classification"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/tem_classification/analysis"
default_metab_table <- file.path(if (.Platform$OS.type == "windows") "C:/projects" else "/mnt/c/projects",
                                 "trec_biobank/metab/analysis/tables/06_trophic_modes_site_fraction.tsv")

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(patchwork)
})
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)
metab_table <- Sys.getenv("TEM_METAB_TROPHIC_TABLE", default_metab_table)
metab_sample_table <- Sys.getenv("TEM_METAB_SAMPLE_TABLE", file.path(dirname(metab_table), "05_trophic_mode_per_sample.tsv"))

## ---- settings ----
metab_fractions <- c("5-40", "10-40", "10-200", "20-200")   # metaB size fractions shown
metab_set       <- "dinoflagellates"                        # or "all eukaryotes"

## metab site name, TEM site code, coordinates (WGS84), bar shift dx / dy (m)
site_coords <- tribble(
  ~site,                       ~tem,  ~lon,     ~lat,     ~dx,    ~dy,
  "Athens",                    "ATH", 23.7275,  37.9838,  0,      0,
  "Barcelona",                 "BAR",  2.1734,  41.3851,  1.5e5,  -4.5e5,
  "Bilbao",                    "BIL", -2.9349,  43.2630,  0,      0,
  "Iceland 2022",              NA,   -21.8277,  64.1283,  0,      0,
  "Kristineberg",              "KRI", 11.4394,  58.2375,  0,      0,
  "Naples",                    "NAP", 14.2681,  40.8522,  0,      0,
  "Porto",                     "POR", -8.6110,  41.1496,  -3.5e5, -6.5e5,
  "Roscoff",                   NA,    -3.9833,  48.7269,  -4e5,   1e5,
  "Tallinn",                   "TAL", 24.7536,  59.4370,  0,      0,
  "Vigo",                      NA,    -8.7207,  42.2406,  -4e5,   1.5e5,
  "Villefranche-sur-Mer 2021", NA,     7.3120,  43.7040,  0,      0
)

tem_levels <- c("chloroplast + food vacuole", "chloroplast only", "no chloroplast")
## aligned with metaB: both chloroplast groups green (= phototroph, light / dark), no chloroplast = heterotroph orange
tem_cols   <- c(`chloroplast + food vacuole` = "#7CC8B0", `chloroplast only` = "#1B9E77", `no chloroplast` = "#D95F02")
metab_levels <- c("phototroph", "photosymbiosis", "heterotroph", "Metazoa", "NA: genus not in W4", "NA: no genus (PR2)")
metab_cols   <- c(phototroph = "#1B9E77", photosymbiosis = "#A6D854", heterotroph = "#D95F02",
                  Metazoa = "#5B7DB1", `NA: genus not in W4` = "grey55", `NA: no genus (PR2)` = "grey82")

map_crs <- "+proj=ortho +lat_0=52 +lon_0=5"
bar_w   <- 0.9e5    # width of one bar (m)
bar_gap <- 0.15e5   # gap between bars; between the TEM bar and the metaB bars 3x
bar_h   <- 4e5      # bar height (100 %)
bar_w_frac   <- 0.7e5    # by-fraction version: narrower bars (more bars per site)
bar_gap_frac <- 0.08e5
bar_h_frac   <- 3.3e5    # lower bars, so that the labels of neighbouring sites do not collide
margin  <- 3.5e5    # map margin around the sites

## ---- data ----
tem_file  <- file.path(opt$table_dir, "chloro_foodvac_site_counts.tsv")
cell_file <- file.path(opt$table_dir, "chloro_foodvac_cells.tsv")
for (f in c(tem_file, cell_file)) if (!file.exists(f)) stop("Missing ", f, ": run 11_chloro_foodvac_sites.R first")
for (f in c(metab_table, metab_sample_table)) if (!file.exists(f))
  stop("Missing ", f, ": run metab/analysis/05 + 06 first (or set TEM_METAB_TROPHIC_TABLE / TEM_METAB_SAMPLE_TABLE)")

tem_site_name <- function(code) {
  out <- site_coords$site[match(code, site_coords$tem)]
  if (any(is.na(out))) stop("TEM site(s) without coordinates: ", paste(unique(code[is.na(out)]), collapse = ", "))
  out
}
frac_label <- function(x) str_replace(x, "to", "-")          # TEM "10to40" -> "10-40"
frac_order <- function(x) {                                  # by lower, then upper bound
  lo <- as.numeric(str_extract(x, "^[\\d.]+")); hi <- as.numeric(str_extract(x, "[\\d.]+$"))
  rank(lo * 1e4 + hi, ties.method = "min")
}

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

## ---- base map ----
sf_use_s2(FALSE)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  suppressMessages(st_crop(., xmin = -35, xmax = 45, ymin = 28, ymax = 72)) %>%
  st_transform(map_crs)
sf_use_s2(TRUE)

sites_all <- site_coords %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(map_crs)
sites_all <- bind_cols(st_drop_geometry(sites_all), as_tibble(st_coordinates(sites_all))) %>%
  mutate(bx = X + dx, by = Y + dy)

## ---- map with bars ----
## bars_long: site, source, bar (id within the site), pair (bars of one pair stand closer together),
##            bar_order, label, n, cat, share
map_with_bars <- function(bars_long, title, subtitle, site_labels, w = bar_w, gap = bar_gap, bar_h = get("bar_h", envir = globalenv())) {
  sites <- sites_all %>% filter(site %in% bars_long$site)
  slots <- bars_long %>%
    distinct(site, bar, pair, source, label, bar_order) %>%
    arrange(site, bar_order) %>%
    group_by(site) %>%
    mutate(step = w + gap + if_else(lag(pair, default = first(pair)) != pair, 2 * gap, 0),
           pos = cumsum(step) - first(step), pos = pos - (min(pos) + max(pos)) / 2) %>%
    ungroup() %>%
    left_join(sites %>% select(site, bx, by), by = "site") %>%
    mutate(xc = bx + pos)
  bars <- bars_long %>%
    left_join(slots %>% select(site, bar, xc, by), by = c("site", "bar")) %>%
    arrange(site, bar, desc(cat)) %>%                         # first category on top
    group_by(site, bar) %>%
    mutate(top = cumsum(share), bottom = lag(top, default = 0)) %>%
    ungroup() %>%
    mutate(xmin = xc - w / 2, xmax = xc + w / 2, ymin = by + bottom * bar_h, ymax = by + top * bar_h)
  tem_frames <- bars %>% filter(source == "TEM") %>% group_by(site, bar, xmin, xmax) %>%
    summarise(ymin = min(ymin), ymax = max(ymax), .groups = "drop")
  groups <- slots %>% group_by(site) %>%
    summarise(xmin = min(xc) - w / 2, xmax = max(xc) + w / 2, .groups = "drop") %>%
    left_join(sites, by = "site") %>%
    mutate(label = site_labels[site])
  moved <- groups %>% filter(dx != 0 | dy != 0)
  xlim <- range(c(sites$X, groups$xmin, groups$xmax)) + c(-1, 1) * margin
  ylim <- c(min(c(sites$Y, sites$by)) - margin, max(sites$by) + bar_h + margin * 0.6)

  p <- ggplot() +
    geom_sf(data = world, fill = "#E4E4E0", colour = "white", linewidth = 0.25) +
    ## line from the real position to the nearest corner of the moved bar group
    geom_segment(data = moved, aes(x = X, y = Y, xend = case_when(dx < 0 ~ xmax, dx > 0 ~ xmin, TRUE ~ bx),
                                   yend = if_else(dy < 0, by + bar_h, by)),
                 colour = "grey35", linewidth = 0.3) +
    geom_segment(data = groups, aes(x = xmin, xend = xmax, y = by, yend = by), colour = "grey35", linewidth = 0.3) +
    geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cat),
              colour = "white", linewidth = 0.15) +
    geom_rect(data = tem_frames, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, colour = "black", linewidth = 0.45) +
    geom_text(data = slots, aes(xc, by - 3.5e4, label = label, fontface = if_else(source == "TEM", "bold", "plain")),
              angle = 90, hjust = 1, size = 1.9, colour = "grey15") +
    geom_point(data = sites, aes(X, Y), shape = 21, fill = "black", colour = "white", size = 1.6, stroke = 0.4) +
    geom_text(data = groups, aes((xmin + xmax) / 2, by + bar_h, label = label), vjust = -0.3, size = 2.6, lineheight = 0.9) +
    scale_fill_manual(values = c(tem_cols, metab_cols), guide = "none") +
    coord_sf(crs = map_crs, xlim = xlim, ylim = ylim, expand = FALSE, datum = NA) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(panel.background = element_rect(fill = "#F3F7FA", colour = NA),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "grey30", size = 8))
  list(plot = p, bars = bars, asp = diff(ylim) / diff(xlim))
}

## two legends (TEM groups / metaB trophic modes) from dummy plots
legend_of <- function(bars, levels, cols, title, nrow) {
  used <- intersect(levels, as.character(unique(bars$cat[bars$share > 0])))
  p <- ggplot(tibble(cat = factor(used, used)), aes(cat, fill = cat)) + geom_bar() +
    scale_fill_manual(values = cols[used], name = title) +
    guides(fill = guide_legend(nrow = nrow, byrow = TRUE)) +
    theme_void(base_size = 11) + theme(legend.position = "bottom", legend.title.position = "top",
                                       legend.key.spacing.x = unit(14, "pt"))
  cowplot::get_plot_component(p, "guide-box-bottom")
}
legends_of <- function(bars) {
  wrap_elements(legend_of(bars, tem_levels, tem_cols, "TEM: dinoflagellate cells", nrow = 3)) |
    wrap_elements(legend_of(bars, metab_levels, metab_cols, "metaB: trophic mode of the reads (Tara W4)", nrow = 2))
}
metab_reads <- if (metab_set == "dinoflagellates") "dinoflagellate" else "18S"
white_bg <- theme(plot.background = element_rect(fill = "white", colour = NA))

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
  site_labels = setNames(site_coords$site, site_coords$site))
p1 <- m1$plot / legends_of(m1$bars) + plot_layout(heights = c(1, 0.08)) & white_bg
save_plot(p1, "trophic_map_combined", 11, 11 * m1$asp + 1.6, opt$plot_dir)
write_tsv(bars1 %>% select(site, source, size_fraction = frac, n, cat, share), file.path(opt$table_dir, "trophic_map_combined.tsv"))

## ---- 2) TEM and metaB per size fraction + sampling time span ----
fmt_span <- function(d) {
  a <- min(d); b <- max(d)
  if (format(a, "%Y%m") == format(b, "%Y%m")) {
    if (a == b) format(a, "%-d %b %Y") else paste0(format(a, "%-d"), "-", format(b, "%-d %b %Y"))
  } else if (format(a, "%Y") == format(b, "%Y")) paste0(format(a, "%-d %b"), " - ", format(b, "%-d %b %Y"))
  else paste0(format(a, "%-d %b %Y"), " - ", format(b, "%-d %b %Y"))
}
fmt_months <- function(d) {
  a <- min(d); b <- max(d)
  if (format(a, "%Y%m") == format(b, "%Y%m")) format(a, "%b %Y")
  else if (format(a, "%Y") == format(b, "%Y")) paste0(format(a, "%b"), "-", format(b, "%b %Y"))
  else paste0(format(a, "%b %Y"), " - ", format(b, "%b %Y"))
}
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
  site_labels = site_labels2, w = bar_w_frac, gap = bar_gap_frac, bar_h = bar_h_frac)

## timeline: one panel per site (sorted by sampling date), TEM and metaB as two lines;
## all panels cover the same number of days (one day = same length in every panel),
## starting 1 day before the first sampling day of the site; panel label = site + month(s) + year
site_order <- dates %>% group_by(site) %>%
  summarise(first = min(date), last = max(date), months = fmt_months(date), .groups = "drop") %>%
  arrange(first) %>%
  mutate(year = format(first, "%Y"),
         months = if_else(str_detect(site, year), str_trim(str_remove(months, year)), months),
         panel = paste0(site, "  ", months))
n_days <- as.numeric(max(site_order$last - site_order$first)) + 2     # longest campaign + 1 day each side
panel_limits <- site_order %>%
  transmute(site, date = first - 1) %>%
  bind_rows(site_order %>% transmute(site, date = first - 1 + n_days))
spans <- dates %>%
  group_by(site, source) %>%
  summarise(start = min(date), end = max(date), n_dates = n_distinct(date), .groups = "drop")
time_data <- function(d) d %>%
  left_join(site_order %>% select(site, panel), by = "site") %>%
  mutate(panel = factor(panel, site_order$panel))
day_labels <- function(x) {                                  # day; month added at the first break / month start
  x <- as.Date(x); ok <- !is.na(x)
  first <- ok & (x == min(x[ok]) | format(x, "%d") %in% c("01", "02"))
  ifelse(first, format(x, "%d\n%b"), format(x, "%d"))
}
src_cols <- c(TEM = "black", metaB = "grey50")
p_time <- ggplot(time_data(spans)) +
  geom_blank(data = time_data(panel_limits) %>% mutate(source = "TEM"), aes(date, source)) +
  geom_segment(aes(x = start, xend = end, y = source, yend = source, colour = source), linewidth = 2.6, alpha = 0.25, lineend = "round") +
  geom_point(data = time_data(dates), aes(date, source, colour = source), size = 1.8) +
  geom_text(aes(x = end, y = source, label = n_dates, colour = source), nudge_y = 0.38, size = 2.1) +
  facet_wrap(~ panel, ncol = 1, scales = "free_x", strip.position = "top") +
  scale_colour_manual(values = src_cols, guide = "none") +
  scale_x_date(breaks = scales::breaks_width("2 days"), labels = day_labels, expand = expansion(0)) +
  scale_y_discrete(limits = c("metaB", "TEM"), expand = expansion(add = 0.7)) +
  labs(x = NULL, y = NULL, title = "Sampling dates",
       subtitle = sprintf(paste0("black = TEM (annotated cells), grey = metaB (shown size fractions)\n",
                                 "number = sampling days; same time scale in all panels (%d days)"), n_days)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "grey97", colour = NA),
        strip.text = element_text(hjust = 0, face = "bold", size = 8.5),
        axis.text.x = element_text(size = 6.5, lineheight = 0.85), panel.spacing.y = unit(3, "pt"),
        plot.title.position = "plot", plot.subtitle = element_text(colour = "grey30", size = 7.5))

p2 <- (m2$plot | p_time) + plot_layout(widths = c(3, 1))
p2 <- p2 / legends_of(m2$bars) + plot_layout(heights = c(1, 0.08)) & white_bg
## map panel width ~ (figure width - axis labels of the timeline) * 3/4
save_plot(p2, "trophic_map_combined_by_fraction", 16, (16 - 1) * 3 / 4 * m2$asp + 2.4, opt$plot_dir)
write_tsv(bars2 %>% select(site, source, size_fraction = frac, n, cat, share),
          file.path(opt$table_dir, "trophic_map_combined_by_fraction.tsv"))
