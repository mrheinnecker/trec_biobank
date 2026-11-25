## Map of Europe with the chloroplast / food vacuole groups of 11_chloro_foodvac_sites.R per site
## one stacked bar (fraction of the three groups) per site, standing on the site's position
## (orthographic projection centred on Europe, as in analysis/metab/ov_map_plot.R)
## Input:  tables/chloro_foodvac_site_counts.tsv of 11_chloro_foodvac_sites.R (run 11 first)
## Output: plots/chloro_foodvac_map

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- NA_character_   # not read here (only the tables of 11)
default_outdir     <- "C:/projects/trec_biobank/analysis/tem/rstudio"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/trec_biobank/analysis/tem"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
})
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)

## ---- settings ----
## site code -> city and coordinates (city centre / marine station, WGS84)
site_coords <- tribble(
  ~site, ~city,          ~lon,     ~lat,
  "ATH", "Athens",       23.7275,  37.9838,
  "BAR", "Barcelona",     2.1734,  41.3851,
  "BIL", "Bilbao",       -2.9349,  43.2630,
  "KRI", "Kristineberg", 11.4394,  58.2375,
  "NAP", "Naples",       14.2681,  40.8522,
  "POR", "Porto",        -8.6110,  41.1496,
  "TAL", "Tallinn",      24.7536,  59.4370
)
## the same groups / colours as 11_chloro_foodvac_sites.R
group_levels <- c("chloroplast + food vacuole", "chloroplast only", "no chloroplast")
group_cols   <- c(`chloroplast + food vacuole` = "#D95F02",
                  `chloroplast only`           = "#1B9E77",
                  `no chloroplast`             = "#7570B3")
map_crs   <- "+proj=ortho +lat_0=50 +lon_0=10"
bar_w     <- 1.7e5   # bar width in projection units (m)
bar_h     <- 4.5e5   # bar height (100 %)
margin    <- 4e5     # map margin around the sites

## ---- data ----
counts_file <- file.path(opt$table_dir, "chloro_foodvac_site_counts.tsv")
if (!file.exists(counts_file)) stop("Missing ", counts_file, ": run 11_chloro_foodvac_sites.R first")
counts <- read_tsv(counts_file, show_col_types = FALSE) %>%
  mutate(troph = factor(troph, group_levels))

unknown <- setdiff(unique(counts$site), site_coords$site)
if (length(unknown) > 0) stop("No coordinates for site(s): ", paste(unknown, collapse = ", "),
                              " -> add them to site_coords")

## site positions in the map projection
sites <- site_coords %>%
  filter(site %in% counts$site) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(map_crs)
sites <- bind_cols(st_drop_geometry(sites), as_tibble(st_coordinates(sites))) %>%
  left_join(distinct(counts, site, n_site), by = "site")

## one stacked bar per site, standing on the site's position
bars <- counts %>%
  left_join(sites %>% select(site, X, Y), by = "site") %>%
  arrange(site, desc(troph)) %>%                       # first group on top
  group_by(site) %>%
  mutate(top = cumsum(frac), bottom = lag(top, default = 0)) %>%
  ungroup() %>%
  mutate(xmin = X - bar_w / 2, xmax = X + bar_w / 2,
         ymin = Y + bottom * bar_h, ymax = Y + top * bar_h,
         label = ifelse(n > 0 & frac >= 0.12, n, ""))

## ---- base map ----
## countries cut to a box around Europe before projecting (the ortho projection cannot show
## the far side of the globe)
sf_use_s2(FALSE)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  suppressMessages(st_crop(., xmin = -30, xmax = 50, ymin = 25, ymax = 75)) %>%
  st_transform(map_crs)
sf_use_s2(TRUE)

xlim <- range(sites$X) + c(-1, 1) * margin
ylim <- c(min(sites$Y) - margin, max(sites$Y) + bar_h + margin * 0.6)

## ---- plot ----
p_map <- ggplot() +
  geom_sf(data = world, fill = "#E4E4E0", colour = "white", linewidth = 0.25) +
  geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = troph),
            colour = "white", linewidth = 0.3) +
  geom_text(data = bars, aes(x = X, y = (ymin + ymax) / 2, label = label), colour = "white", size = 2.6) +
  geom_point(data = sites, aes(X, Y), shape = 21, fill = "black", colour = "white", size = 1.8, stroke = 0.4) +
  geom_text(data = sites, aes(X, Y, label = sprintf("%s (%s)\nn = %d", city, site, n_site)),
            vjust = 1.25, size = 2.8, lineheight = 0.9) +
  scale_fill_manual(values = group_cols, name = NULL, drop = FALSE) +
  coord_sf(crs = map_crs, xlim = xlim, ylim = ylim, expand = FALSE, datum = NA) +
  labs(title = "Chloroplast / food vacuole groups across sites",
       subtitle = "dinoflagellates; bar = fraction of the cells of a site, numbers = cells",
       x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(panel.background = element_rect(fill = "#F3F7FA", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "bottom", plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey30", size = 9),
        plot.margin = margin(8, 8, 8, 8))

asp <- diff(ylim) / diff(xlim)
save_plot(p_map, "chloro_foodvac_map", 9, 9 * asp + 1, opt$plot_dir)
write_tsv(sites, file.path(opt$table_dir, "chloro_foodvac_map_sites.tsv"))
