# Map of Europe with one stacked bar per site x size fraction, for 06_ and 07_.
# Same map as trec_biobank/analysis/tem/12_chloro_foodvac_map.R (orthographic projection,
# Natural Earth countries). The bars of a site stand side by side on the site's position
# (size fractions small to large); sites that are too close to each other are moved by
# `dx` / `dy` (map units, m) with a line to the real position.
# Source after common.R.

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
map_crs <- "+proj=ortho +lat_0=52 +lon_0=5"
bar_w   <- 0.9e5    # width of one bar (m)
bar_gap <- 0.15e5   # gap between the bars of a site
bar_h   <- 4e5      # bar height (100 %)
margin  <- 3.5e5    # map margin around the sites

# ---- helpers -----------------------------------------------------------------
check_site_coords <- function(sites_in_data) {
  unknown <- setdiff(unique(sites_in_data), site_coords$site)
  if (length(unknown) > 0) stop("No coordinates for site(s): ", paste(unknown, collapse = ", "), " -> add them to site_coords in map_common.R")
}

# numeric bounds of "10-40" etc., to order the size fractions
add_fraction_bounds <- function(d) {
  d %>% mutate(lower = as.numeric(str_extract(size_fraction, "^[\\d.]+")),
               upper = as.numeric(str_extract(size_fraction, "[\\d.]+$")))
}

## site positions in the map projection (+ moved bar position)
site_positions <- function(sites_in_data) {
  s <- site_coords %>%
    filter(site %in% sites_in_data) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform(map_crs)
  bind_cols(st_drop_geometry(s), as_tibble(st_coordinates(s))) %>%
    mutate(bx = X + dx, by = Y + dy)
}

base_map <- function() {
  sf_use_s2(FALSE)
  on.exit(sf_use_s2(TRUE))
  ne_countries(scale = "medium", returnclass = "sf") %>%
    st_make_valid() %>%
    suppressMessages(st_crop(., xmin = -35, xmax = 45, ymin = 28, ymax = 72)) %>%
    st_transform(map_crs)
}

# ---- plot --------------------------------------------------------------------
# d: one row per site x size_fraction x fill level, with columns site, size_fraction,
#    n_samples, share and the column named in `fill` (a factor; first level on top)
bar_map <- function(d, fill, colours, legend_title, title,
                    subtitle = "one bar per size fraction (um; number of samples): mean share of the reads over the samples",
                    world = base_map(), legend_nrow = NULL) {
  d <- add_fraction_bounds(d) %>% rename(fill_level = all_of(fill))
  sites <- site_positions(unique(d$site))
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
    arrange(site, size_fraction, desc(fill_level)) %>%        # first level on top
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
    geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_level),
              colour = "white", linewidth = 0.15) +
    ## size fraction (+ number of samples) under every bar
    geom_text(data = slots, aes(xc, by - 3.5e4, label = sprintf("%s (%d)", size_fraction, n_samples)),
              angle = 90, hjust = 1.05, size = 1.9, colour = "grey20") +
    geom_point(data = sites, aes(X, Y), shape = 21, fill = "black", colour = "white", size = 1.6, stroke = 0.4) +
    geom_text(data = groups, aes((xmin + xmax) / 2, by + bar_h, label = site), vjust = -0.5, size = 2.6) +
    scale_fill_manual(values = colours, name = legend_title, drop = FALSE,
                      guide = guide_legend(nrow = legend_nrow)) +
    coord_sf(crs = map_crs, xlim = xlim, ylim = ylim, expand = FALSE, datum = NA) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(panel.background = element_rect(fill = "#F3F7FA", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.position = "bottom", plot.title.position = "plot",
          plot.subtitle = element_text(colour = "grey30", size = 9),
          plot.margin = margin(8, 8, 8, 8))
}

save_map <- function(p, name, width = 11) {
  b <- ggplot_build(p)$layout$panel_params[[1]]
  asp <- diff(b$y_range) / diff(b$x_range)
  save_plot(p, name, width = width, height = width * asp + 1.2)
}
