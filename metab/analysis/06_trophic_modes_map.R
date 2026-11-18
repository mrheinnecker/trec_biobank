# 06: Map of Europe with the trophic modes of 05_trophic_modes.R per site and size fraction
#
# One stacked bar per site x size fraction = mean share of the trophic modes over the samples
# of that site and size fraction (every sample weighted equally). The bars of a site stand
# side by side on the site's position (size fractions small to large); sites that are too
# close to each other are moved by `dx` / `dy` (map units, m) with a line to the real position.
# Same map as tem_classification/analysis/12_chloro_foodvac_map.R (orthographic projection,
# Natural Earth countries).
# Input:  tables/05_trophic_mode_per_sample.tsv of 05_trophic_modes.R (run 05 first)
# Output: plots/06_trophic_modes_map_dinoflagellates, plots/06_trophic_modes_map_eukaryotes,
#         tables/06_trophic_modes_site_fraction.tsv

## hard-coded for running line by line (RStudio); with Rscript the script's own folder is used
default_script_dir <- "C:/repos/tem_classification/metab/analysis"
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
})

# ---- settings ----------------------------------------------------------------
# site (as in sample_mapping.tsv) -> coordinates (WGS84); dx / dy move the bars (m)
site_coords <- tribble(
  ~site,                       ~short, ~lon,     ~lat,     ~dx,    ~dy,
  "Athens",                    "ATH",  23.7275,  37.9838,  0,      0,
  "Barcelona",                 "BAR",   2.1734,  41.3851,  2e5,    -5e5,
  "Bilbao",                    "BIL",  -2.9349,  43.2630,  0,      0,
  "Iceland 2022",              "ICE", -21.8277,  64.1283,  0,      0,
  "Kristineberg",              "KRI",  11.4394,  58.2375,  0,      0,
  "Naples",                    "NAP",  14.2681,  40.8522,  0,      0,
  "Porto",                     "POR",  -8.6110,  41.1496,  -3e5,   -5.5e5,
  "Roscoff",                   "ROS",  -3.9833,  48.7269,  0,      0,
  "Tallinn",                   "TAL",  24.7536,  59.4370,  0,      0,
  "Vigo",                      "VIG",  -8.7207,  42.2406,  -4e5,   0,
  "Villefranche-sur-Mer 2021", "VIL",   7.3120,  43.7040,  0,      0
)
# from 05_trophic_modes.R
mode_levels <- c("phototroph", "photosymbiosis", "heterotroph", "Metazoa",
                 "NA: genus not in W4", "NA: no genus (PR2)")
mode_colours <- c(phototroph = "#1B9E77", photosymbiosis = "#A6D854", heterotroph = "#D95F02",
                  Metazoa = "#5B7DB1", `NA: genus not in W4` = "grey55", `NA: no genus (PR2)` = "grey82")
map_crs <- "+proj=ortho +lat_0=52 +lon_0=5"
bar_w   <- 0.9e5    # width of one bar (m)
bar_gap <- 0.15e5   # gap between the bars of a site
bar_h   <- 4e5      # bar height (100 %)
margin  <- 3.5e5    # map margin around the sites

# ---- data --------------------------------------------------------------------
shares_f <- file.path(table_dir, "05_trophic_mode_per_sample.tsv")
if (!file.exists(shares_f)) stop("Missing ", shares_f, ": run 05_trophic_modes.R first")
shares <- read_tsv(shares_f, show_col_types = FALSE) %>%
  filter(sample_type == "environmental", !is.na(site), !is.na(size_fraction))

unknown <- setdiff(unique(shares$site), site_coords$site)
if (length(unknown) > 0) stop("No coordinates for site(s): ", paste(unknown, collapse = ", "), " -> add them to site_coords")

## mean share per site x size fraction x mode; size fractions ordered by lower, then upper bound
site_fraction <- shares %>%
  group_by(set, site, size_fraction, mode) %>%
  summarise(share = mean(share), n_samples = n_distinct(run_accession), .groups = "drop") %>%
  mutate(mode = factor(mode, mode_levels),
         lower = as.numeric(str_extract(size_fraction, "^[\\d.]+")),
         upper = as.numeric(str_extract(size_fraction, "[\\d.]+$")))
write_tsv(site_fraction %>% select(-lower, -upper), file.path(table_dir, "06_trophic_modes_site_fraction.tsv"))

## site positions in the map projection (+ moved bar position)
sites <- site_coords %>%
  filter(site %in% shares$site) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(map_crs)
sites <- bind_cols(st_drop_geometry(sites), as_tibble(st_coordinates(sites))) %>%
  mutate(bx = X + dx, by = Y + dy)

# ---- base map ----------------------------------------------------------------
sf_use_s2(FALSE)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  suppressMessages(st_crop(., xmin = -35, xmax = 45, ymin = 28, ymax = 72)) %>%
  st_transform(map_crs)
sf_use_s2(TRUE)

# ---- plot --------------------------------------------------------------------
map_plot <- function(set_name, title) {
  d <- site_fraction %>% filter(set == set_name)
  ## one slot per site x size fraction, centred on the (moved) site position
  slots <- d %>%
    distinct(site, size_fraction, lower, upper, n_samples) %>%
    arrange(site, lower, upper) %>%
    group_by(site) %>%
    mutate(i = row_number(), k = n(), offset = (i - (k + 1) / 2) * (bar_w + bar_gap)) %>%
    ungroup() %>%
    left_join(sites %>% select(site, bx, by), by = "site") %>%
    mutate(xc = bx + offset)
  bars <- d %>%
    left_join(slots %>% select(site, size_fraction, xc, by), by = c("site", "size_fraction")) %>%
    arrange(site, size_fraction, desc(mode)) %>%            # first mode on top
    group_by(site, size_fraction) %>%
    mutate(top = cumsum(share), bottom = lag(top, default = 0)) %>%
    ungroup() %>%
    mutate(xmin = xc - bar_w / 2, xmax = xc + bar_w / 2, ymin = by + bottom * bar_h, ymax = by + top * bar_h)
  groups <- slots %>% group_by(site) %>%
    summarise(xmin = min(xc) - bar_w / 2, xmax = max(xc) + bar_w / 2, .groups = "drop") %>%
    left_join(sites, by = "site")
  moved <- groups %>% filter(dx != 0 | dy != 0)

  xlim <- range(c(sites$X, groups$xmin, groups$xmax)) + c(-1, 1) * margin
  ylim <- c(min(c(sites$Y, sites$by)) - margin, max(sites$by) + bar_h + margin * 0.6)

  ggplot() +
    geom_sf(data = world, fill = "#E4E4E0", colour = "white", linewidth = 0.25) +
    geom_segment(data = moved, aes(x = X, y = Y, xend = bx, yend = if_else(dy < 0, by + bar_h, by)),
                 colour = "grey35", linewidth = 0.3) +   # to the top of the bars if they were moved down
    geom_segment(data = groups, aes(x = xmin, xend = xmax, y = by, yend = by), colour = "grey35", linewidth = 0.3) +
    geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = mode),
              colour = "white", linewidth = 0.15) +
    ## size fraction (+ number of samples) under every bar
    geom_text(data = slots, aes(xc, by - 3.5e4, label = sprintf("%s (%d)", size_fraction, n_samples)),
              angle = 90, hjust = 1.05, size = 1.9, colour = "grey20") +
    geom_point(data = sites, aes(X, Y), shape = 21, fill = "black", colour = "white", size = 1.6, stroke = 0.4) +
    geom_text(data = groups, aes((xmin + xmax) / 2, by + bar_h, label = site), vjust = -0.5, size = 2.6) +
    scale_fill_manual(values = mode_colours, name = "trophic mode (W4)", drop = FALSE) +
    coord_sf(crs = map_crs, xlim = xlim, ylim = ylim, expand = FALSE, datum = NA) +
    labs(title = title,
         subtitle = "one bar per size fraction (um; number of samples): mean share of the reads over the samples",
         x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(panel.background = element_rect(fill = "#F3F7FA", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.position = "bottom", plot.title.position = "plot",
          plot.subtitle = element_text(colour = "grey30", size = 9),
          plot.margin = margin(8, 8, 8, 8))
}

save_map <- function(p, name) {
  b <- ggplot_build(p)$layout$panel_params[[1]]
  asp <- diff(b$y_range) / diff(b$x_range)
  save_plot(p, name, width = 11, height = 11 * asp + 1.2)
}
save_map(map_plot("dinoflagellates", "Trophic modes of the dinoflagellate reads (PR2 -> Tara W4)"),
         "06_trophic_modes_map_dinoflagellates")
save_map(map_plot("all eukaryotes", "Trophic modes of the 18S reads (PR2 -> Tara W4)"),
         "06_trophic_modes_map_eukaryotes")
cat("plots written to", plot_dir, "\n")
