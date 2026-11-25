## Shared map code of 13_trophic_map_combined.R and 14_tem_feature_maps.R:
## site coordinates (with bar shifts), base map of Europe, stacked bars on the map, legends,
## date formatting and the sampling-date timeline. Needs tidyverse, sf, rnaturalearth, patchwork.

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


map_crs <- "+proj=ortho +lat_0=52 +lon_0=5"
bar_w   <- 0.9e5    # width of one bar (m)
bar_gap <- 0.15e5   # gap between bars; between the TEM bar and the metaB bars 3x
bar_h   <- 4e5      # bar height (100 %)
bar_w_frac   <- 0.7e5    # by-fraction version: narrower bars (more bars per site)
bar_gap_frac <- 0.08e5
bar_h_frac   <- 3.3e5    # lower bars, so that the labels of neighbouring sites do not collide
margin  <- 3.5e5    # map margin around the sites

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
##            bar_order, label, n, cat (factor), share; bars with source "TEM" get a black outline
##            (and a bold label) -> use another source name if no outline is wanted
map_with_bars <- function(bars_long, title, subtitle, site_labels, fill_values, w = bar_w, gap = bar_gap,
                          bar_h = get("bar_h", envir = globalenv())) {
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
    summarise(ymin = min(ymin), ymax = max(ymax), .groups = "drop") %>% suppressWarnings()   # no TEM bars -> empty
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
    scale_fill_manual(values = fill_values, guide = "none") +
    coord_sf(crs = map_crs, xlim = xlim, ylim = ylim, expand = FALSE, datum = NA) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(panel.background = element_rect(fill = "#F3F7FA", colour = NA),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "grey30", size = 8))
  list(plot = p, bars = bars, asp = diff(ylim) / diff(xlim))
}

## legend of the used categories (levels -> colours) from a dummy plot
legend_of <- function(bars, levels, cols, title, nrow) {
  used <- intersect(levels, as.character(unique(bars$cat[bars$share > 0])))
  p <- ggplot(tibble(cat = factor(used, used)), aes(cat, fill = cat)) + geom_bar() +
    scale_fill_manual(values = cols[used], name = title) +
    guides(fill = guide_legend(nrow = nrow, byrow = TRUE)) +
    theme_void(base_size = 11) + theme(legend.position = "bottom", legend.title.position = "top",
                                       legend.key.spacing.x = unit(14, "pt"))
  cowplot::get_plot_component(p, "guide-box-bottom")
}
white_bg <- theme(plot.background = element_rect(fill = "white", colour = NA))

## ---- dates ----
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

## timeline: one panel per site (sorted by sampling date), one line per source (TEM / metaB);
## all panels cover the same number of days (one day = same length in every panel),
## starting 1 day before the first sampling day of the site; panel label = site + month(s) + year
timeline_plot <- function(dates, subtitle_sources = "black = TEM (annotated cells), grey = metaB (shown size fractions)") {
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
ggplot(time_data(spans)) +
  geom_blank(data = time_data(panel_limits) %>% mutate(source = "TEM"), aes(date, source)) +
  geom_segment(aes(x = start, xend = end, y = source, yend = source, colour = source), linewidth = 2.6, alpha = 0.25, lineend = "round") +
  geom_point(data = time_data(dates), aes(date, source, colour = source), size = 1.8) +
  geom_text(aes(x = end, y = source, label = n_dates, colour = source), nudge_y = 0.38, size = 2.1) +
  facet_wrap(~ panel, ncol = 1, scales = "free_x", strip.position = "top") +
  scale_colour_manual(values = src_cols, guide = "none") +
  scale_x_date(breaks = scales::breaks_width("2 days"), labels = day_labels, expand = expansion(0)) +
  scale_y_discrete(limits = intersect(c("metaB", "TEM"), unique(dates$source)), expand = expansion(add = 0.7)) +
  labs(x = NULL, y = NULL, title = "Sampling dates",
       subtitle = sprintf("%s\nnumber = sampling days; same time scale in all panels (%d days)", subtitle_sources, n_days)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = "grey97", colour = NA),
        strip.text = element_text(hjust = 0, face = "bold", size = 8.5),
        axis.text.x = element_text(size = 6.5, lineheight = 0.85), panel.spacing.y = unit(3, "pt"),
        plot.title.position = "plot", plot.subtitle = element_text(colour = "grey30", size = 7.5))
}
