## Biobank overview, step 2: figures (A circle, B HPF flow, C radial flow, D radial chain).
## Input: a data snapshot of biobank_overview.R (<fig1_root>/data/<timestamp>/, default the newest).
## Every run writes into a new folder <fig1_root>/runs/<YYYYmmdd_HHMMSS>/ with the plots, a copy of
## this script (to compare the settings of two runs) and run_info.txt (data snapshot, git commit).

library(tidyverse)
library(circlize)
library(ggalluvial)

## ---- run folder ----
fig1_root <- if (.Platform$OS.type == "windows") "C:/projects/trec_biobank/fig1" else "/mnt/c/projects/trec_biobank/fig1"

## data snapshot to use (folder name in <fig1_root>/data); NA = the newest
## (env FIG1_DATA_SNAPSHOT, set by launcher.sh, overrides this)
data_snapshot <- NA
if (Sys.getenv("FIG1_DATA_SNAPSHOT") != "") data_snapshot <- Sys.getenv("FIG1_DATA_SNAPSHOT")
ov_dir <- if (is.na(data_snapshot)) {
  sort(list.dirs(file.path(fig1_root, "data"), recursive=FALSE), decreasing=TRUE)[1]
} else {
  file.path(fig1_root, "data", data_snapshot)
}
if (is.na(ov_dir) || !file.exists(file.path(ov_dir, "event_coverage.tsv"))) {
  stop("No data snapshot with event_coverage.tsv found (run biobank_overview.R first): ", ov_dir)
}

## optional note appended to the run folder name, e.g. "new_colours" -> 20261118_170000_new_colours
## (env FIG1_RUN_NOTE overrides this)
run_note <- ""
if (Sys.getenv("FIG1_RUN_NOTE") != "") run_note <- Sys.getenv("FIG1_RUN_NOTE")
plot_dir <- file.path(fig1_root, "runs", paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), if (nzchar(run_note)) paste0("_", run_note)))
dir.create(plot_dir, recursive=TRUE)
message("Data: ", ov_dir, "\nOutput: ", plot_dir)

## run info + a copy of this script (Rscript, source() or the open file in RStudio)
script_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error=function(e) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
    if (is.na(f)) NA_character_ else normalizePath(f)
  }
)
if (is.na(script_file) && requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  script_file <- rstudioapi::getSourceEditorContext()$path
}
git_commit <- tryCatch(
  system2("git", c("-C", dirname(script_file), "rev-parse", "--short", "HEAD"), stdout=TRUE, stderr=FALSE),
  error=function(e) NA_character_, warning=function(w) NA_character_
)
writeLines(c(
  paste("time:", format(Sys.time())),
  paste("data snapshot:", ov_dir),
  paste("script:", script_file),
  paste("git commit:", paste(git_commit, collapse=""))
), file.path(plot_dir, "run_info.txt"))
if (!is.na(script_file) && file.exists(script_file)) file.copy(script_file, plot_dir)

## a PDF that is open in a viewer cannot be overwritten on Windows: warn instead of stopping
ggsave_pdf <- function(filename, plot, ...) {
  tryCatch(
    ggsave(filename, plot, ...),
    error=function(e) warning(basename(filename), " not written (open in a viewer?)")
  )
}


## =================================================================================================
## ---- COLOURS AND LEGEND LABELS: edit here -------------------------------------------------------
## =================================================================================================
## One row per element: key = name in the data (do not change), color = fill colour,
## label = text in the legend (NA = left out of the legend). The legend lists the entries in row
## order. Two entries with the same colour share one legend entry.

## sites: site ring and site-coloured flows of all figures; the row order is the site order
site_palette <- tribble(
  ~key,           ~color,    ~label,
  "Roscoff",      "#332288", NA,
  "Tallinn",      "#4477AA", NA,
  "Kristineberg", "#66CCEE", NA,
  "Bilbao",       "#228833", NA,
  "Vigo",         "#44AA99", NA,
  "Porto",        "#AA9933", NA,
  "Barcelona",    "#EE7733", NA,
  "Naples",       "#CC6677", NA,
  "Athens",       "#AA3377", NA
)

## D (chain figure): blocks of rings 2-5 (the flows outside the site ring take the colour of the
## block they come from, lightened by outer_flow_lighten further down)
chain_palette <- tribble(
  ~key,       ~color,    ~label,
  ## ring 2: samples
  "HPF",      "#6EC6CA", "cryo-storage",
  "PF",       "#6EC6CA", "cryo-storage",
  "MetaB",    "#5AAE61", "ENA",
  "Chem fix", "#4477AA", "4°C-storage",
  ## ring 3: resin
  "Epon",     "#FDAE6B", "RT-storage",
  "Lowicryl", "#FDAE6B", "RT-storage",
  ## ring 4: imaging of the sample
  "HiTT",     "#35978F", "BioImageArchive",
  "resin LM", "#35978F", "BioImageArchive",
  "cryo LM",  "#35978F", "BioImageArchive",
  ## ring 5: cells
  "TEM",      "#35978F", "BioImageArchive",
  "SEM",      "#35978F", "BioImageArchive",
  "ExM",      "#35978F", "BioImageArchive"
)
## legend of D: title, position ("right", "bottom", ..., "none" = no legend), extra figure width (in)
## for a legend on the right, and whether the sites are listed too (their label column above)
chain_legend_title <- NULL
chain_legend_position <- "right"
chain_legend_width <- 1.6
chain_legend_sites <- FALSE

## E (radial chain split by site): blocks of rings 2-5, same keys as chain_palette;
## ring = legend block the entry is listed in (titles in radial_chain_ring_titles);
## repository = where the samples / data end up, shown as symbol(s) next to the legend entry
## (several separated by ";", NA = no symbol; symbols in repository_palette)
radial_chain_palette <- tribble(
  ~ring, ~key,       ~color,    ~label,                      ~repository,
  2,     "HPF",      "#2166AC", "HPF",                       "BioSamples",
  2,     "PF",       "#4393C3", "plunge frozen (PF + WM)",   "BioSamples",   # only used with e_split_pf <- FALSE
  2,     "cryo PF",  "#4393C3", "plunge frozen (PF)",        "BioSamples",
  2,     "cryo WM",  "#92C5DE", "waffle method (WM)",        "BioSamples",
  2,     "MetaB",    "#5AAE61", "metabarcoding",             "BioSamples",
  2,     "Chem fix", "#9970AB", "chemically fixed",          "BioSamples",
  3,     "Epon",     "#B2182B", "Epon",                      "BioSamples",
  3,     "Lowicryl", "#EF8A62", "Lowicryl",                  "BioSamples",
  4,     "V4V5 rRNA sequencing", "#1B7837", "V4V5 rRNA sequencing", "ENA",   # ring = e_metab_sequencing_ring
  4,     "HiTT",     "#D94801", "HiTT (X-ray)",              "BioImageArchive",
  4,     "resin LM", "#FD8D3C", "resin LM",                  "BioImageArchive",
  4,     "cryo LM",  "#FDAE6B", "cryo LM",                   "BioImageArchive",
  5,     "TEM",      "#7F2704", "TEM",                       "BioImageArchive",
  5,     "SEM",      "#8C6D31", "SEM",                       "BioImageArchive",
  5,     "ExM",      "#35978F", "ExM",                       "BioImageArchive"
)
## repositories: symbol (filled point shape: 21 circle, 22 square, 23 diamond, 24 triangle up,
## 25 triangle down), colour and the text in the repository key; the key lists them in row order
repository_palette <- tribble(
  ~repository,       ~shape, ~color,    ~label,
  "BioSamples",      21,     "#4D4D4D", "BioSamples",
  "ENA",             24,     "#E69F00", "ENA",
  "BioImageArchive", 22,     "#0072B2", "BioImageArchive"
)
repository_legend_title <- "repository"
## legend of E: one block per ring, from inside to outside, with these titles
radial_chain_ring_titles <- c(
  `2`="raw sample",
  `3`="resin blocks",
  `4`="overview scan/dataset",
  `5`="images"
)
## legend: "right" (drawn next to the circle, with the repository symbols) or "none"; whether the
## sites are listed too (as a first block with this title)
radial_chain_legend_position <- "right"
radial_chain_legend_sites <- FALSE
radial_chain_sites_title <- "site"

## C (radial flow): sample ring and image ring
sample_type_palette <- tribble(
  ~key,       ~color,    ~label,
  "HPF",      "#2166AC", "HPF",
  "cryo PF",  "#4393C3", "cryo PF",
  "cryo WM",  "#92C5DE", "cryo WM",
  "MetaB",    "#5AAE61", "MetaB",
  "Chem fix", "#9970AB", "Chem fix"
)
imaging_palette <- tribble(
  ~key,          ~color,    ~label,
  "TEM",         "#7F2704", "TEM",
  "HiTT",        "#D94801", "HiTT",
  "resin LM",    "#FD8D3C", "resin LM",
  "cryo LM PF",  "#FDAE6B", "cryo LM PF",
  "cryo LM WM",  "#FDD0A2", "cryo LM WM",
  "cryo LM HPF", "#FEE6CE", "cryo LM HPF",
  "SEM",         "#8C6D31", "SEM",
  "ExM",         "#35978F", "ExM"
)
event_only_color <- "grey75"
event_only_label <- "no sample record"

## A (circle): colour ramps for 0 / 1 / 10 / 100 / 1000 samples or images per event
sample_ramp <- c("grey94", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
image_ramp <- c("grey94", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704")
## =================================================================================================

palette_colors <- function(pal) deframe(select(pal, key, color))
## legend entries (colour -> label) of the rows with a label, in row order, one per colour
palette_legend <- function(pal) {
  pal %>%
    filter(!is.na(label)) %>%
    distinct(color, .keep_all=TRUE) %>%
    with(setNames(label, color))
}

site_colors <- palette_colors(site_palette)

coverage <-read_tsv(file.path(ov_dir, "event_coverage.tsv"), show_col_types=FALSE)
hpf_events <- read_tsv(file.path(ov_dir, "hpf_event_map.tsv"), col_types=cols(.default=col_character()))
imaging_map <- read_tsv(file.path(ov_dir, "imaging_map.tsv"), col_types=cols(.default=col_character()))


## ---- A: circular overview, one column per sampling event ----------------------------------------
## sectors = sites, events in time order within a site;
## inner rings = number of samples / images per event and modality (log scale)

sample_rings <- c(
  HPF="n_hpf",
  `cryo PF`="n_cryo_pf",
  `cryo WM`="n_cryo_wm",
  MetaB="n_metab",
  `Chem fix`="n_chemfix"
)
image_rings <- c(
  TEM="n_img_tem",
  HiTT="n_img_hitt",
  `resin LM`="n_img_resin_lm",
  `cryo LM PF`="n_img_cryo_lm_pf",
  `cryo LM WM`="n_img_cryo_lm_wm",
  `cryo LM HPF`="n_img_cryo_lm_hpf",
  SEM="n_img_sem",
  ExM="n_img_exm"
)

circle_data <- coverage %>%
  mutate(
    site=factor(site, levels=names(site_colors)),
    across(any_of(c(sample_rings, image_rings)), ~coalesce(as.numeric(.x), 0))
  ) %>%
  arrange(site, date, event_time)

## modalities without any image yet (missing column) count as 0
for (col in setdiff(c(sample_rings, image_rings), names(circle_data))) {
  circle_data[[col]] <- 0
}

ring_matrix <- function(rings) {
  m <- circle_data %>%
    select(all_of(unname(rings))) %>%
    as.matrix()
  colnames(m) <- names(rings)
  log10(m + 1)
}

sample_matrix <- ring_matrix(sample_rings)
image_matrix <- ring_matrix(image_rings)

## log10(n + 1): 0 = none, 0.3 = 1, 1 = 9, 2 = 99, 3 = 999
count_breaks <- c(0, log10(2), 1, 2, 3)
sample_col_fun <- colorRamp2(count_breaks, sample_ramp)
image_col_fun <- colorRamp2(count_breaks, image_ramp)

campaign_labels <- circle_data %>%
  group_by(site) %>%
  summarise(
    campaign=if_else(
      format(min(date), "%b %Y") == format(max(date), "%b %Y"),
      format(min(date), "%b %Y"),
      paste(format(min(date), "%b"), "-", format(max(date), "%b %Y"))
    )
  ) %>%
  deframe()

split <- droplevels(circle_data$site)
n_sectors <- nlevels(split)

## ring names in the gap after the last sector (circos.heatmap draws the first column outermost)
add_ring_labels <- function(labels) {
  circos.track(
    track.index=get.current.track.index(),
    panel.fun=function(x, y) {
      if (CELL_META$sector.numeric.index == n_sectors) {
        n <- length(labels)
        circos.text(
          rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
          seq_len(n) - 0.5,
          rev(labels),
          cex=0.45,
          adj=c(0, 0.5),
          facing="inside",
          niceFacing=TRUE
        )
      }
    },
    bg.border=NA
  )
}

draw_circle <- function() {
  circos.clear()
  circos.par(
    start.degree=90,
    gap.after=c(rep(2, n_sectors - 1), 16),
    cell.padding=c(0, 0, 0, 0),
    track.margin=c(0.004, 0.004)
  )

  ## site ring with site name and sampling campaign
  circos.heatmap(
    matrix(as.character(circle_data$site), ncol=1),
    split=split,
    col=site_colors,
    track.height=0.04,
    cluster=FALSE
  )
  circos.track(
    track.index=get.current.track.index(),
    panel.fun=function(x, y) {
      circos.text(
        CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(5),
        CELL_META$sector.index,
        facing="bending.inside", niceFacing=TRUE, cex=0.8, font=2
      )
      circos.text(
        CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(1.5),
        campaign_labels[[CELL_META$sector.index]],
        facing="bending.inside", niceFacing=TRUE, cex=0.5
      )
    },
    bg.border=NA
  )

  circos.heatmap(sample_matrix, split=split, col=sample_col_fun, track.height=0.2, cluster=FALSE, bg.border="grey60", bg.lwd=0.5)
  add_ring_labels(colnames(sample_matrix))

  circos.heatmap(image_matrix, split=split, col=image_col_fun, track.height=0.24, cluster=FALSE, bg.border="grey60", bg.lwd=0.5)
  add_ring_labels(colnames(image_matrix))

  ## totals in the middle
  n_samples <- sum(circle_data[unname(sample_rings)])
  n_images <- sum(circle_data[unname(image_rings)])
  text(
    0, 0,
    paste0(
      nrow(circle_data), " sampling events\n",
      format(n_samples, big.mark=","), " samples\n",
      format(n_images, big.mark=","), " images"
    ),
    cex=0.8
  )

  legend_counts <- c(0, 1, 10, 100, 1000)
  legend(
    "bottomleft", inset=c(0, 0), bty="n", cex=0.6, title="samples per event",
    legend=c("none", "1", "10", "100", "1000"),
    fill=sample_col_fun(log10(legend_counts + 1)), border="grey60"
  )
  legend(
    "bottomright", inset=c(0, 0), bty="n", cex=0.6, title="images per event",
    legend=c("none", "1", "10", "100", "1000"),
    fill=image_col_fun(log10(legend_counts + 1)), border="grey60"
  )
  circos.clear()
}

png(file.path(plot_dir, "biobank_overview_circle.png"), width=3000, height=3000, res=400)
draw_circle()
dev.off()

pdf_ok <- tryCatch(
  {
    pdf(file.path(plot_dir, "biobank_overview_circle.pdf"), width=7.5, height=7.5)
    TRUE
  },
  error=function(e) {
    warning("biobank_overview_circle.pdf not written (open in a viewer?)")
    FALSE
  }
)
if (pdf_ok) {
  draw_circle()
  dev.off()
}


## ---- B: HPF samples: site -> processing -> imaging -----------------------------------------------
## block_name = embedded in resin, cryo_name = kept in cryo

hpf_imaging <- imaging_map %>%
  filter(sample_layer == "HPF") %>%
  distinct(sample_barcode, imaging) %>%
  group_by(barcode=sample_barcode) %>%
  summarise(imaging=paste(sort(imaging), collapse=" + "))

hpf_flow <- hpf_events %>%
  mutate(
    site=factor(site, levels=names(site_colors)),
    processing=case_when(
      block_name == "discarded" ~ "discarded",
      block_name == "block?" ~ "resin (unclear)",
      !is.na(block_name) ~ "resin block",
      !is.na(cryo_name) ~ "cryo",
      TRUE ~ "unknown"
    )
  ) %>%
  left_join(hpf_imaging, by="barcode") %>%
  mutate(
    ## cryo HPF: imaged with cryo LM according to the sheet, but not (yet) in the collection table
    imaging=case_when(
      !is.na(imaging) ~ imaging,
      processing == "cryo" & imaged == "TRUE" ~ "cryo LM HPF",
      TRUE ~ "not imaged"
    )
  ) %>%
  count(site, processing, imaging)

processing_levels <- c("resin block", "resin (unclear)", "cryo", "discarded", "unknown")
imaging_levels <- c(
  hpf_flow %>% filter(imaging != "not imaged") %>% count(imaging, wt=n, sort=TRUE) %>% pull(imaging),
  "not imaged"
)

## small imaging strata: listed next to the imaging axis, starting at the top of the first small one
## (ggalluvial stacks the strata in level order from the top)
small_label_min <- 40
imaging_counts <- hpf_flow %>%
  count(imaging, wt=n) %>%
  mutate(imaging=factor(imaging, levels=imaging_levels)) %>%
  arrange(imaging) %>%
  mutate(top=sum(n) - lag(cumsum(n), default=0))
small_imaging <- imaging_counts %>% filter(n < small_label_min)

flow_plot <- hpf_flow %>%
  mutate(
    processing=factor(processing, levels=processing_levels),
    imaging=factor(imaging, levels=imaging_levels)
  ) %>%
  ggplot(aes(axis1=site, axis2=processing, axis3=imaging, y=n)) +
  geom_alluvium(aes(fill=site), alpha=0.75, width=0.35) +
  geom_stratum(width=0.35, fill="grey97", color="grey40", linewidth=0.3) +
  ## large strata labelled inside
  geom_text(
    stat="stratum",
    aes(label=if_else(after_stat(count) >= small_label_min, paste0(after_stat(stratum), "\n", after_stat(count)), "")),
    size=2.4, lineheight=0.9
  ) +
  annotate(
    "text",
    x=3.22,
    y=if (nrow(small_imaging) > 0) max(small_imaging$top) else 0,
    label=paste0(small_imaging$imaging, " (", small_imaging$n, ")", collapse="\n"),
    hjust=0, vjust=1, size=2.2, lineheight=0.95
  ) +
  scale_x_discrete(limits=c("Site", "Processing", "Imaging"), expand=expansion(add=c(0.25, 0.7))) +
  scale_fill_manual(values=site_colors) +
  labs(x=NULL, y=NULL, title=paste0(format(sum(hpf_flow$n), big.mark=","), " HPF samples")) +
  theme_minimal(base_size=9) +
  theme(
    panel.grid=element_blank(),
    axis.text.y=element_blank(),
    legend.position="none"
  )

ggsave(file.path(plot_dir, "biobank_overview_hpf_flow.png"), flow_plot, width=7, height=5, dpi=300, bg="white")
ggsave_pdf(file.path(plot_dir, "biobank_overview_hpf_flow.pdf"), flow_plot, width=7, height=5)


## ---- C: radial flow: sites -> samples -> images -----------------------------------------------------
## inner ring: sites sized by sampling events, middle ring: samples per site and type,
## outer ring: images per site and imaging method; flows bundle all samples / images of one path.
## Each ring fills the circle in its own unit (events, samples, images), so a flow can be narrower
## at one end: at the sample end it is as wide as the number of imaged samples (the not imaged rest
## of a sample block has no outgoing flow). Images linked to an event only (no sample record, e.g. SEM)
## go through a grey "no sample record" block in the sample ring, sized in sample units.

## colours: sample_type_palette / imaging_palette in the colour section at the top
sample_type_colors <- palette_colors(sample_type_palette)
imaging_colors <- palette_colors(imaging_palette)
event_only <- "event only"

sample_maps <- c(
  HPF="hpf_event_map.tsv",
  `cryo PF`="cryo_pf_event_map.tsv",
  `cryo WM`="cryo_wm_event_map.tsv",
  MetaB="metab_event_map.tsv",
  `Chem fix`="chemfix_event_map.tsv"
)

site_levels <- names(site_colors)

ring1_nodes <- coverage %>%
  filter(site %in% site_levels) %>%
  count(site, name="size") %>%
  mutate(site=factor(site, levels=site_levels)) %>%
  arrange(site)

samples_per_type <- imap_dfr(sample_maps, function(file, sample_type) {
  read_tsv(file.path(ov_dir, file), col_types=cols(.default=col_character())) %>%
    ## as in the circle: samples with a sampling event
    filter(site %in% site_levels, !is.na(event_id)) %>%
    count(site, name="size") %>%
    mutate(sample_type=sample_type)
})

image_paths <- imaging_map %>%
  filter(!is.na(event_id)) %>%
  inner_join(coverage %>% select(event_id, site), by="event_id") %>%
  filter(site %in% site_levels) %>%
  mutate(sample_type=recode(coalesce(sample_layer, event_only), `Chem fixation`="Chem fix")) %>%
  group_by(site, sample_type, imaging) %>%
  summarise(n_images=n(), n_samples=n_distinct(sample_barcode, na.rm=TRUE), .groups="drop")

## the gap for event-only images is sized in sample units: images * (all samples / all images)
samples_per_image <- sum(samples_per_type$size) / sum(image_paths$n_images)

ring2_nodes <- bind_rows(
  samples_per_type,
  image_paths %>%
    filter(sample_type == event_only) %>%
    group_by(site, sample_type) %>%
    summarise(size=sum(n_images) * samples_per_image, .groups="drop")
) %>%
  mutate(
    site=factor(site, levels=site_levels),
    sample_type=factor(sample_type, levels=c(names(sample_type_colors), event_only))
  ) %>%
  arrange(site, sample_type)

ring3_nodes <- image_paths %>%
  group_by(site, imaging) %>%
  summarise(size=sum(n_images), .groups="drop") %>%
  mutate(
    site=factor(site, levels=site_levels),
    imaging=factor(imaging, levels=names(imaging_colors))
  ) %>%
  arrange(site, imaging)

## angles in degrees, clockwise from the top; gaps between sites and between blocks of one site
site_gap <- 3
node_gap <- 0.6

layout_ring <- function(nodes) {
  n_nodes <- nrow(nodes)
  n_sites <- n_distinct(nodes$site)
  span <- 360 - n_sites * site_gap - (n_nodes - n_sites) * node_gap
  nodes %>%
    mutate(
      width=size / sum(size) * span,
      new_site=site != lag(site, default=first(site)) | row_number() == 1,
      gap=if_else(new_site, site_gap, node_gap),
      ## half a site gap before the first site, so the circle is symmetric around the top
      gap=if_else(row_number() == 1, site_gap / 2, gap),
      start=cumsum(gap) + lag(cumsum(width), default=0),
      end=start + width
    ) %>%
    select(-new_site, -gap)
}

ring1 <- layout_ring(ring1_nodes)
ring2 <- layout_ring(ring2_nodes)
ring3 <- layout_ring(ring3_nodes)

## pieces of a node for its flows (in node order), widths in node units, filled from the start
split_node <- function(start, end, size, widths) {
  scale <- (end - start) / max(size, sum(widths))
  ends <- start + cumsum(widths) * scale
  tibble(from=c(start, head(ends, -1)), to=ends)
}

## radii
r1 <- c(0.24, 0.32)
r2 <- c(0.52, 0.58)
r3 <- c(0.80, 0.86)

to_xy <- function(angle, r) {
  tibble(x=r * sin(angle * pi / 180), y=r * cos(angle * pi / 180))
}
arc_xy <- function(a1, a2, r, n=40) {
  to_xy(seq(a1, a2, length.out=n), r)
}
block_poly <- function(a1, a2, radii) {
  bind_rows(arc_xy(a1, a2, radii[1]), arc_xy(a2, a1, radii[2]))
}
## flow from [a1, a2] at radius ra to [b1, b2] at radius rb, smooth in angle
## shortest signed angle from a to b
ang_diff <- function(a, b) ((b - a + 180) %% 360) - 180

flow_poly <- function(a1, a2, b1, b2, ra, rb, n=50) {
  ## the shorter way around the circle
  shift <- (a1 + a2) / 2 + ang_diff((a1 + a2) / 2, (b1 + b2) / 2) - (b1 + b2) / 2
  b1 <- b1 + shift
  b2 <- b2 + shift
  t <- seq(0, 1, length.out=n)
  s <- 3 * t^2 - 2 * t^3
  bind_rows(
    to_xy(a1 + (b1 - a1) * s, ra + (rb - ra) * t),
    arc_xy(b1, b2, rb),
    to_xy(a2 + (b2 - a2) * rev(s), ra + (rb - ra) * rev(t)),
    arc_xy(a2, a1, ra)
  )
}

## flows site -> sample block: the site arc is split by the sizes of its sample blocks
flows12 <- ring2 %>%
  group_by(site) %>%
  group_modify(function(blocks, key) {
    site_node <- ring1 %>% filter(site == key$site)
    pieces <- split_node(site_node$start, site_node$end, sum(blocks$size), blocks$size)
    tibble(sample_type=blocks$sample_type, a1=pieces$from, a2=pieces$to, b1=blocks$start, b2=blocks$end)
  }) %>%
  ungroup()

## flows sample block -> image block: width at the sample end = imaged samples (event only: its gap),
## at the image end = images; pieces ordered by imaging method / sample type
flows23 <- image_paths %>%
  mutate(
    site=factor(site, levels=site_levels),
    sample_type=factor(sample_type, levels=levels(ring2$sample_type)),
    imaging=factor(imaging, levels=names(imaging_colors)),
    sample_width=if_else(sample_type == event_only, n_images * samples_per_image, as.numeric(n_samples))
  ) %>%
  arrange(site, sample_type, imaging) %>%
  group_by(site, sample_type) %>%
  group_modify(function(paths, key) {
    node <- ring2 %>% filter(site == key$site, sample_type == key$sample_type)
    pieces <- split_node(node$start, node$end, node$size, paths$sample_width)
    mutate(paths, a1=pieces$from, a2=pieces$to)
  }) %>%
  ungroup() %>%
  arrange(site, imaging, sample_type) %>%
  group_by(site, imaging) %>%
  group_modify(function(paths, key) {
    node <- ring3 %>% filter(site == key$site, imaging == key$imaging)
    pieces <- split_node(node$start, node$end, node$size, paths$n_images)
    mutate(paths, b1=pieces$from, b2=pieces$to)
  }) %>%
  ungroup()

## flows ending in the event-only gap continue straight through the sample ring
flow_polygons <- bind_rows(
  flows12 %>%
    mutate(id=paste("f12", row_number()), fill=site_colors[as.character(site)], alpha=0.55) %>%
    rowwise() %>%
    reframe(id, fill, alpha, flow_poly(a1, a2, b1, b2, r1[2], r2[1])),
  flows23 %>%
    mutate(id=paste("f23", row_number()), fill=site_colors[as.character(site)], alpha=0.55) %>%
    rowwise() %>%
    reframe(id, fill, alpha, flow_poly(a1, a2, b1, b2, r2[2], r3[1])),

)

block_polygons <- bind_rows(
  ring1 %>%
    mutate(id=paste("r1", row_number()), fill=site_colors[as.character(site)]) %>%
    rowwise() %>%
    reframe(id, fill, block_poly(start, end, r1)),
  ring2 %>%
    mutate(
      id=paste("r2", row_number()),
      fill=if_else(sample_type == event_only, event_only_color, sample_type_colors[as.character(sample_type)])
    ) %>%
    rowwise() %>%
    reframe(id, fill, block_poly(start, end, r2)),
  ring3 %>%
    mutate(id=paste("r3", row_number()), fill=imaging_colors[as.character(imaging)]) %>%
    rowwise() %>%
    reframe(id, fill, block_poly(start, end, r3))
)

## text along the circle (tangential) or pointing outwards (radial), flipped to stay readable
tangential_angle <- function(mid) if_else(mid > 90 & mid < 270, 180 - mid, -mid)
radial_angle <- function(mid) if_else(mid < 180, 90 - mid, 270 - mid)

site_labels <- ring1 %>%
  mutate(
    mid=(start + end) / 2,
    x=to_xy(mid, mean(r1))$x,
    y=to_xy(mid, mean(r1))$y,
    angle=tangential_angle(mid),
    ## short sites: name only
    label=if_else(end - start >= 25, paste0(site, "\n", size, " events"), as.character(site))
  )

image_labels <- ring3 %>%
  filter(end - start >= 2.5) %>%
  mutate(
    mid=(start + end) / 2,
    x=to_xy(mid, r3[2] + 0.015)$x,
    y=to_xy(mid, r3[2] + 0.015)$y,
    angle=radial_angle(mid),
    hjust=if_else(mid < 180, 0, 1),
    label=paste0(imaging, " (", size, ")")
  )

sample_labels <- ring2 %>%
  filter(sample_type != event_only, end - start >= 4) %>%
  mutate(
    mid=(start + end) / 2,
    x=to_xy(mid, mean(r2))$x,
    y=to_xy(mid, mean(r2))$y,
    angle=tangential_angle(mid),
    label=size
  )

## legend: samples, the event-only block, images (labels from the palettes at the top)
radial_legend <- c(
  palette_legend(sample_type_palette),
  if (!is.na(event_only_label)) setNames(event_only_label, event_only_color),
  palette_legend(imaging_palette)
)
radial_legend <- radial_legend[!duplicated(names(radial_legend))]

radial_plot <- ggplot() +
  geom_polygon(data=flow_polygons, aes(x, y, group=id, fill=fill), alpha=0.5, color=NA) +
  geom_polygon(data=block_polygons, aes(x, y, group=id, fill=fill), color="white", linewidth=0.15) +
  geom_text(data=site_labels, aes(x, y, label=label, angle=angle), color="white", size=2, lineheight=0.85, fontface="bold") +
  geom_text(data=sample_labels, aes(x, y, label=label, angle=angle), color="white", size=1.7) +
  geom_text(data=image_labels, aes(x, y, label=label, angle=angle, hjust=hjust), size=1.7) +
  annotate(
    "text", 0, 0, size=2.6, lineheight=0.95,
    label=paste0(
      sum(ring1$size), " sampling events\n",
      format(sum(samples_per_type$size), big.mark=","), " samples\n",
      format(sum(ring3$size), big.mark=","), " images"
    )
  ) +
  scale_fill_identity(
    guide="legend",
    breaks=names(radial_legend),
    labels=unname(radial_legend),
    name=NULL
  ) +
  coord_equal(xlim=c(-1.05, 1.05), ylim=c(-1.05, 1.05)) +
  theme_void(base_size=8) +
  theme(
    legend.position="right",
    legend.key.size=unit(3, "mm"),
    plot.background=element_rect(fill="white", color=NA)
  )

ggsave(file.path(plot_dir, "biobank_overview_radial.png"), radial_plot, width=8, height=6.6, dpi=400, bg="white")
ggsave_pdf(file.path(plot_dir, "biobank_overview_radial.pdf"), radial_plot, width=8, height=6.6)


## helpers for D
layout_blocks <- function(sizes, deg_per_unit=NULL) {
  sizes <- sizes[sizes > 0]
  if (is.null(deg_per_unit)) {
    deg_per_unit <- (360 - length(sizes) * block_gap) / sum(sizes)
  }
  widths <- sizes * deg_per_unit
  start <- block_gap / 2 + cumsum(c(0, head(widths + block_gap, -1)))
  tibble(block=names(sizes), size=unname(sizes), start=start, end=start + widths)
}

best_rotation <- function(src_center, dst_center, weight) {
  offsets <- 0:359
  cost <- map_dbl(offsets, ~sum(weight * abs(ang_diff(src_center, dst_center + .x))))
  offsets[which.min(cost)]
}

rotate_blocks <- function(blocks, offset) mutate(blocks, start=start + offset, end=end + offset)
block_centers <- function(blocks) transmute(blocks, block, center=(start + end) / 2)
site_center <- function(site) {
  idx <- match(as.character(site), as.character(ring1$site))
  (ring1$start[idx] + ring1$end[idx]) / 2
}

block_centers <- function(blocks) transmute(blocks, block, center=(start + end) / 2)
site_center <- function(site) {
  idx <- match(as.character(site), as.character(ring1$site))
  (ring1$start[idx] + ring1$end[idx]) / 2
}

site_center <- function(site) {
  idx <- match(as.character(site), as.character(ring1$site))
  (ring1$start[idx] + ring1$end[idx]) / 2
}


## ---- D: radial chain with merged rings: sites -> samples -> resin blocks -> imaging -> cells -------
## ring 1: sites (sized by events)
## ring 2: sample types (samples)
## ring 3: HPF samples embedded in resin, Epon / Lowicryl from the block name (_epo / _low) (samples)
## ring 4: imaging of the physical sample: HiTT, resin LM, cryo LM (images)
## ring 5: "cells": TEM, SEM and ExM (usually one cell per image)
## Rings 2-5 are not split by site, the site is carried by the flow colour. Paths that skip a ring
## (PF -> cryo LM, chem fix -> SEM, resin block -> TEM without HiTT) run through an invisible
## pass-through block of that ring, drawn as one continuous flow.
## TEM comes from HiTT when the same HPF sample (block) has HiTT images.
## Rings 2 / 3 share one scale (degrees per sample), rings 4 / 5 one scale (degrees per image).

## TRUE: the flows outside the site ring keep the site colour (one flow per site);
## FALSE: one flow per path, light grey (outer_flow_color)
color_flows_by_site <- FALSE

## colour of the flows outside the site ring when they are not coloured by site
outer_flow_color <- "grey85"
## TRUE: the flows outside the site ring take the colour of the block they come from (e.g. HPF -> Epon
## in the HPF blue), mixed with white by outer_flow_lighten (0 = block colour, 1 = white);
## FALSE: outer_flow_color
color_flows_by_source <- TRUE
outer_flow_lighten <- 0.4
lighten <- function(color, amount) {
  rgb_mix <- col2rgb(color) * (1 - amount) + 255 * amount
  rgb(rgb_mix[1, ], rgb_mix[2, ], rgb_mix[3, ], maxColorValue=255)
}

## blocks without _epo / _low in the name (block?, a few without suffix) count as Epon
resin_types <- c("Epon", "Lowicryl")
resin_of <- function(x) {
  if_else(str_detect(x, regex("_low", ignore_case=TRUE)), "Lowicryl", "Epon")
}
pass_prefix <- "pass:"
is_pass <- function(x) str_starts(x, fixed(pass_prefix))

## colours: chain_palette in the colour section at the top
chain_block_colors <- palette_colors(chain_palette)
ring2_levels <- c("HPF", "PF", "MetaB", "Chem fix")

merge_sample_type <- function(x) {
  ## plunge frozen: cryo PF and WM together
  recode(x, `cryo PF`="PF", `cryo WM`="PF", `Chem fixation`="Chem fix")
}

chain_samples <- samples_per_type %>%
  mutate(sample_type=merge_sample_type(sample_type)) %>%
  count(site, sample_type, wt=size, name="size")

chain_images <- imaging_map %>%
  filter(!is.na(event_id), !is.na(sample_layer)) %>%
  inner_join(coverage %>% select(event_id, site), by="event_id") %>%
  filter(site %in% site_levels) %>%
  mutate(
    sample_type=merge_sample_type(sample_layer),
    imaging=recode(imaging, `cryo LM PF`="cryo LM", `cryo LM WM`="cryo LM", `cryo LM HPF`="cryo LM"),
    sample_key=coalesce(sample_key, image_id)
  )

## the ExM images come from chemically fixed FA2 aliquots, which are not in the chem fixation sheet:
## they are added to the chem fix samples
chain_samples <- bind_rows(
  chain_samples,
  chain_images %>%
    filter(imaging == "ExM") %>%
    distinct(site, sample_key) %>%
    count(site, name="size") %>%
    mutate(sample_type="Chem fix")
) %>%
  count(site, sample_type, wt=size, name="size")

## resin per HPF sample: block name in the HPF sheet, else the block / image name of its images
## resin blocks not in the HPF sheet, linked to the event only (Vigo log; written by biobank_overview.R)
blocks_extra_file <- file.path(ov_dir, "resin_blocks_extra.tsv")
blocks_extra <- if (file.exists(blocks_extra_file)) {
  read_tsv(blocks_extra_file, col_types=cols(.default=col_character())) %>% filter(!is.na(event_id))
} else {
  tibble(site=character(), event_id=character(), resin=character(), sample_key=character())
}

key_resin <- bind_rows(
  blocks_extra %>% transmute(sample_key, resin),
  hpf_events %>%
    filter(!is.na(block_name), block_name != "discarded") %>%
    transmute(sample_key=barcode, resin=resin_of(block_name)),
  chain_images %>%
    filter(sample_type == "HPF", imaging %in% c("TEM", "HiTT", "resin LM")) %>%
    transmute(sample_key, resin=resin_of(coalesce(block_name, image_id)))
) %>%
  distinct(sample_key, .keep_all=TRUE)

chain_images <- chain_images %>%
  left_join(key_resin, by="sample_key")

## HPF samples in resin: blocks in the HPF sheet and blocks only known from TEM / HiTT / resin LM images
embedded <- bind_rows(
  hpf_events %>%
    filter(!is.na(event_id), site %in% site_levels, !is.na(block_name), block_name != "discarded") %>%
    transmute(site, sample_key=barcode),
  chain_images %>%
    filter(sample_type == "HPF", imaging %in% c("TEM", "HiTT", "resin LM")) %>%
    distinct(site, sample_key),
  blocks_extra %>% filter(site %in% site_levels) %>% transmute(site, sample_key)
) %>%
  distinct(site, sample_key) %>%
  left_join(key_resin, by="sample_key") %>%
  mutate(resin=coalesce(resin, "Epon"))

keys_with <- function(img) {
  chain_images %>% filter(imaging == img) %>% distinct(sample_key) %>% pull()
}
hitt_keys <- keys_with("HiTT")
tem_keys <- keys_with("TEM")

images_of <- function(img, keys=NULL) {
  out <- chain_images %>% filter(imaging == img)
  if (!is.null(keys)) out <- out %>% filter(sample_key %in% keys)
  out
}
## per site (and resin): samples (distinct keys) and images; "{resin}" in from / to is the resin type
per_site <- function(img_rows, from, to, w_from, w_to) {
  img_rows %>%
    group_by(site, resin) %>%
    summarise(n_keys=n_distinct(sample_key), n_img=n(), .groups="drop") %>%
    transmute(
      site,
      from=str_replace(from, fixed("{resin}"), coalesce(resin, "")),
      to=str_replace(to, fixed("{resin}"), coalesce(resin, "")),
      w_from={{ w_from }},
      w_to={{ w_to }}
    )
}
pass_tem <- paste(pass_prefix, "TEM {resin}")

## segments between ring k and k + 1, widths in the units of both rings
segments <- list(
  `1`=chain_samples %>% transmute(site, from=site, to=sample_type, w_from=size, w_to=size),
  `2`=bind_rows(
    embedded %>% count(site, resin, name="n") %>% transmute(site, from="HPF", to=resin, w_from=n, w_to=n),
    per_site(images_of("cryo LM") %>% filter(sample_type == "PF"), "PF", paste(pass_prefix, "PF"), n_keys, n_keys),
    ## HPF carriers kept in cryo and screened with cryo LM
    per_site(images_of("cryo LM") %>% filter(sample_type == "HPF"), "HPF", paste(pass_prefix, "HPF cryo"), n_keys, n_keys),
    per_site(images_of("SEM"), "Chem fix", paste(pass_prefix, "Chem fix"), n_keys, n_keys),
    per_site(images_of("ExM"), "Chem fix", paste(pass_prefix, "Chem fix ExM"), n_keys, n_keys)
  ),
  `3`=bind_rows(
    per_site(images_of("HiTT"), "{resin}", "HiTT", n_keys, n_img),
    per_site(images_of("resin LM"), "{resin}", "resin LM", n_keys, n_img),
    per_site(images_of("TEM", setdiff(tem_keys, hitt_keys)), "{resin}", pass_tem, n_keys, n_img),
    per_site(images_of("cryo LM") %>% filter(sample_type == "PF"), paste(pass_prefix, "PF"), "cryo LM", n_keys, n_img),
    per_site(images_of("cryo LM") %>% filter(sample_type == "HPF"), paste(pass_prefix, "HPF cryo"), "cryo LM", n_keys, n_img),
    per_site(images_of("SEM"), paste(pass_prefix, "Chem fix"), paste(pass_prefix, "SEM"), n_keys, n_img),
    per_site(images_of("ExM"), paste(pass_prefix, "Chem fix ExM"), paste(pass_prefix, "ExM"), n_keys, n_img)
  ),
  `4`=bind_rows(
    images_of("HiTT", tem_keys) %>%
      count(site, name="n_hitt") %>%
      inner_join(images_of("TEM", hitt_keys) %>% count(site, name="n_tem"), by="site") %>%
      transmute(site, from="HiTT", to="TEM", w_from=n_hitt, w_to=n_tem),
    per_site(images_of("TEM", setdiff(tem_keys, hitt_keys)), pass_tem, "TEM", n_img, n_img),
    per_site(images_of("SEM"), paste(pass_prefix, "SEM"), "SEM", n_img, n_img),
    per_site(images_of("ExM"), paste(pass_prefix, "ExM"), "ExM", n_img, n_img)
  )
) %>%
  imap(~mutate(.x, seg=as.integer(.y), w_from=as.numeric(w_from), w_to=as.numeric(w_to))) %>%
  map(~filter(.x, w_from > 0, w_to > 0))

## per-site segments, kept for E (radial chain split by site)
segments_site <- segments

if (!color_flows_by_site) {
  segments[-1] <- map(segments[-1], function(seg) {
    seg %>%
      group_by(seg, from, to) %>%
      summarise(w_from=sum(w_from), w_to=sum(w_to), .groups="drop") %>%
      mutate(site="all")
  })
}

## ring radii
## rings 2-5 with equal gaps (0.06)
chain_radii <- list(c(0.15, 0.21), c(0.33, 0.38), c(0.44, 0.49), c(0.55, 0.60), c(0.66, 0.71))

block_gap <- 2

## blocks of a ring from the flows ending there; sizes in the ring unit
ring_sizes <- function(seg_in) {
  seg_in %>% count(block=to, wt=w_to, name="size")
}

## pieces of a block for its flows, ordered by the angle of the other end (relative to the block centre)
## and by site; fill = the flows fill the block (site ring: events vs samples)
allocate_side <- function(flows, blocks, node_col, other_center_col, prefix, fill=FALSE) {
  flows %>%
    left_join(blocks %>% select(block, b_start=start, b_end=end, b_size=size), by=setNames("block", node_col)) %>%
    mutate(
      rel=ang_diff((b_start + b_end) / 2, .data[[other_center_col]]),
      site_order=as.integer(factor(site, levels=site_levels)),
      ## fixed order of incoming flows (target side), flows put at the clockwise end (source side)
      fixed_rank=if (prefix == "b") {
        map2_dbl(.data[[node_col]], from, ~coalesce(match(.y, incoming_order_fixed[[.x]] %||% character()), 0L) * 1)
      } else {
        map2_dbl(.data[[node_col]], to, ~coalesce(match(.y, outgoing_order_last[[.x]] %||% character()), 0L) * 1)
      }
    ) %>%
    arrange(.data[[node_col]], fixed_rank, rel, site_order) %>%
    group_by(across(all_of(node_col))) %>%
    group_modify(function(f, key) {
      width_col <- if (prefix == "a") "w_from" else "w_to"
      size <- if (fill) sum(f[[width_col]]) else f$b_size[1]
      pieces <- split_node(f$b_start[1], f$b_end[1], size, f[[width_col]])
      ## outgoing_order_last: moved to the very end of the block (the pieces only fill the imaged part)
      at_end <- prefix == "a" & f$fixed_rank > 0
      if (any(at_end)) {
        shift <- f$b_end[1] - max(pieces$to[at_end])
        pieces$from[at_end] <- pieces$from[at_end] + shift
        pieces$to[at_end] <- pieces$to[at_end] + shift
      }
      f[[paste0(prefix, "1")]] <- pieces$from
      f[[paste0(prefix, "2")]] <- pieces$to
      f
    }) %>%
    ungroup() %>%
    select(-b_start, -b_end, -b_size, -rel, -site_order, -fixed_rank)
}

## blocks placed over the weighted centre of their sources, then pushed apart where they overlap
place_blocks <- function(sizes, source_centers, deg_per_unit) {
  blocks <- sizes %>%
    left_join(source_centers, by="block") %>%
    mutate(width=size * deg_per_unit) %>%
    arrange(center)
  for (iteration in 1:200) {
    n <- nrow(blocks)
    if (n < 2) break
    moved <- FALSE
    for (i in seq_len(n)) {
      j <- if (i == n) 1 else i + 1
      need <- (blocks$width[i] + blocks$width[j]) / 2 + block_gap
      ## clockwise distance to the next block
      have <- (blocks$center[j] - blocks$center[i]) %% 360
      if (n == 2 && i == 2) break
      if (have < need) {
        push <- (need - have) / 2
        blocks$center[i] <- blocks$center[i] - push
        blocks$center[j] <- blocks$center[j] + push
        moved <- TRUE
      }
    }
    if (!moved) break
  }
  blocks %>% transmute(block, size, start=center - width / 2, end=center + width / 2)
}

## mean angle of pieces (weights), robust to crossing 0 degrees
circular_center <- function(angles, weights) {
  angles[1] + weighted.mean(ang_diff(angles[1], angles), weights)
}
source_centers <- function(seg) {
  seg %>%
    mutate(mid=(a1 + a2) / 2) %>%
    group_by(block=to) %>%
    summarise(center=circular_center(mid, w_to), .groups="drop")
}

## ring 1: sites as in C
chain_rings <- list(ring1 %>% transmute(block=as.character(site), size, start, end))

## ring 2: sample types around the full circle, rotated onto the sites
## ---- block order optimisation ----
## the order of the blocks in rings 2, 4 and 5 is chosen so that the flows travel the least angle around
## the circle: all orders are tried, each with its best rotation; cost = mean angle travelled per unit of
## flow width, summed over the segments. Ring 1 (sites) keeps its geographic order, ring 3 sits over its
## sources anyway.
permutations <- function(x) {
  if (length(x) <= 1) return(list(x))
  do.call(c, lapply(seq_along(x), function(i) lapply(permutations(x[-i]), function(p) c(x[i], p))))
}
travel <- function(src, dst, w) sum(w * abs(ang_diff(src, dst))) / sum(w)
## best rotation, vectorised over all offsets
best_rotation <- function(src_center, dst_center, weight) {
  offsets <- 0:359
  cost <- vapply(offsets, function(o) sum(weight * abs(ang_diff(src_center, dst_center + o))), numeric(1))
  offsets[which.min(cost)]
}
## blocks in a given order, rotated onto their sources
ordered_ring <- function(sizes, order, deg, src, block, w) {
  sizes <- sizes[match(order, sizes$block), ]
  blocks <- layout_blocks(setNames(sizes$size, sizes$block), deg)
  center <- block_centers(blocks)$center[match(block, blocks$block)]
  rotate_blocks(blocks, best_rotation(src, center, w))
}
ring_cost <- function(blocks, src, block, w) {
  travel(src, block_centers(blocks)$center[match(block, blocks$block)], w)
}

## fixed clockwise order of ring 4 (NULL: optimised); HiTT next to the TEM pass-through above Epon
## fixed clockwise order of ring 5 (NULL: optimised): ExM and SEM swapped against the optimum
ring5_order_fixed <- c("ExM", "SEM", "TEM")
## ring 5 squeezed into a part of the circle (degrees clockwise from 12:00; NULL: full circle, rotated
## onto its sources): the cell imaging in the right half, ExM at 12:00 to TEM at 6:00
ring5_span <- c(0, 180)
ring4_order_fixed <- c("resin LM", "cryo LM", paste(pass_prefix, "SEM"), paste(pass_prefix, "TEM Epon"), "HiTT")

## extra rotation of a ring after the optimisation, in degrees (negative = counter-clockwise);
## the rings outside follow, as their rotation is optimised onto the rotated ring
ring_extra_rotation <- c(`4`=-90, `5`=-20)
## fixed clockwise order of ring 2 (NULL: optimised); Chem fix and MetaB swapped against the optimum
ring2_order_fixed <- c("MetaB", "HPF", "PF", "Chem fix")
## ring 2 rotated so that this block starts at the given angle (NULL: rotated onto the sites);
## HPF from 3:00 on, rings 3 and 4 follow
ring2_anchor <- c(HPF=90)
## manual positions: start angle of a block in degrees clockwise from 12:00 (3:00 = 90, 6:00 = 180,
## 9:00 = 270); the block keeps its width, the flows follow. Applied last, so it overrides the order,
## anchor and rotations above for that block. Invisible pass-through blocks ("pass: ...") can be set
## too, but only reserve space (flows are drawn from their first to their last visible block).
## The current positions are printed on every run.
block_start_fixed <- c(
  ## ring 2: samples
  MetaB=45.4-90,
  HPF=0,
  PF=315.7-90,
  `Chem fix`=11.4-90,
  ## ring 3: resin blocks
  Epon=10,
  Lowicryl=90,
  ## ring 4: imaging
  `resin LM`=95,
  `cryo LM`=200,
  HiTT=10,
  ## ring 5: cells
  ExM=270,
  SEM=23.5-90,
  TEM=10.2
)
## extra rotation of single blocks after their ring is placed, in degrees (positive = clockwise)
block_extra_rotation <- c(`cryo LM`=105, HiTT=20, `resin LM`=-10)
## fixed clockwise order of the flows entering a block (by the block they come from), other flows after
incoming_order_fixed <- list(TEM=c("HiTT", paste(pass_prefix, "TEM Epon")))
## flows leaving a block that are put at its clockwise end (by the block they go to), e.g. the cryo HPF
## flow at the end of HPF next to PF, close to cryo LM
outgoing_order_last <- list(
  HPF=c(paste(pass_prefix, "HPF cryo")),
  ## Chem fix: ExM first, SEM at the end
  `Chem fix`=c(paste(pass_prefix, "Chem fix"))
)

ring2_sizes <- ring_sizes(segments[[1]])
deg_per_sample <- (360 - nrow(ring2_sizes) * block_gap) / sum(ring2_sizes$size)
seg1_src <- segments[[1]] %>% mutate(src=site_center(site))

## ring 2: the first block is fixed (the order is circular), all orders of the others
ring2_candidates <- map(permutations(ring2_sizes$block[-1]), ~c(ring2_sizes$block[1], .x))
ring2_costs <- map_dbl(ring2_candidates, function(order) {
  blocks <- ordered_ring(ring2_sizes, order, deg_per_sample, seg1_src$src, seg1_src$to, seg1_src$w_to)
  ring_cost(blocks, seg1_src$src, seg1_src$to, seg1_src$w_to)
})
ring2_order <- ring2_candidates[[which.min(ring2_costs)]]
if (!is.null(ring2_order_fixed)) {
  fixed <- intersect(ring2_order_fixed, ring2_sizes$block)
  ring2_order <- c(fixed, setdiff(ring2_sizes$block, fixed))
}
message("ring 2 order: ", paste(ring2_order, collapse=", "))
ring2_blocks <- ordered_ring(ring2_sizes, ring2_order, deg_per_sample, seg1_src$src, seg1_src$to, seg1_src$w_to)
if (!is.null(ring2_anchor)) {
  anchor_start <- ring2_blocks$start[ring2_blocks$block == names(ring2_anchor)]
  ring2_blocks <- rotate_blocks(ring2_blocks, ring2_anchor[[1]] - anchor_start)
}
## blocks with a manual start angle (block_start_fixed) moved there
apply_block_start_fixed <- function(blocks) {
  fixed <- block_start_fixed[intersect(names(block_start_fixed), blocks$block)]
  if (length(fixed) == 0) return(blocks)
  blocks %>%
    mutate(
      new_start=coalesce(unname(fixed[block]), start),
      end=end + new_start - start,
      start=new_start
    ) %>%
    select(-new_start)
}
chain_rings[[2]] <- apply_block_start_fixed(ring2_blocks)

## segment 1: site side filled (events vs samples)
segments[[1]] <- segments[[1]] %>%
  mutate(site_mid=site_center(site)) %>%
  left_join(block_centers(ring2_blocks) %>% rename(to_center=center), by=c("to"="block")) %>%
  allocate_side(chain_rings[[1]], "from", "to_center", "a", fill=TRUE) %>%
  allocate_side(chain_rings[[2]], "to", "site_mid", "b")

## rings 3-5: source side first, then the blocks over their sources, then the target side
deg_per_image <- NULL
for (k in 2:4) {
  seg <- segments[[k]]
  ## preliminary target centres: the centre of the source block
  src_block_centers <- block_centers(chain_rings[[k]])
  seg <- seg %>%
    left_join(src_block_centers %>% rename(from_center=center), by=c("from"="block")) %>%
    mutate(to_center_guess=from_center) %>%
    allocate_side(chain_rings[[k]], "from", "to_center_guess", "a")

  sizes <- ring_sizes(seg)
  if (k == 2) {
    deg <- deg_per_sample
  } else {
    if (is.null(deg_per_image)) {
      ## rings 4 and 5 in images: the larger of the two fills the circle
      total4 <- sum(ring_sizes(segments[[3]])$size)
      total5 <- sum(ring_sizes(segments[[4]])$size)
      deg_per_image <- (360 - 6 * block_gap) / max(total4, total5)
    }
    deg <- deg_per_image
  }
  if (k == 2) {
    ## ring 3 (subsets of the ring 2 blocks): each block over its sources
    chain_rings[[k + 1]] <- place_blocks(sizes, source_centers(seg), deg)
  } else {
    ## rings 4 and 5 are almost full: blocks in the optimised order, rotated as a whole
    if (k == 3) {
      ## ring 4 and ring 5 together: all orders of ring 4, for each all orders of ring 5
      src3 <- seg %>% transmute(block=to, w=w_to, src=(a1 + a2) / 2)
      seg4 <- segments[[4]] %>% group_by(from, to) %>% summarise(w=sum(w_to), .groups="drop")
      sizes5 <- ring_sizes(segments[[4]])
      best_cost <- Inf
      order4_candidates <- if (is.null(ring4_order_fixed)) {
        map(permutations(sizes$block[-1]), ~c(sizes$block[1], .x))
      } else {
        ## blocks without flows are left out, blocks missing in the fixed order are appended
        fixed <- intersect(ring4_order_fixed, sizes$block)
        list(c(fixed, setdiff(sizes$block, fixed)))
      }
      for (order4 in order4_candidates) {
        blocks4 <- ordered_ring(sizes, order4, deg, src3$src, src3$block, src3$w)
        cost4 <- ring_cost(blocks4, src3$src, src3$block, src3$w)
        src4 <- block_centers(blocks4)$center[match(seg4$from, blocks4$block)]
        order5_candidates <- if (is.null(ring5_order_fixed)) {
          permutations(sizes5$block)
        } else {
          fixed5 <- intersect(ring5_order_fixed, sizes5$block)
          list(c(fixed5, setdiff(sizes5$block, fixed5)))
        }
        for (order5 in order5_candidates) {
          blocks5 <- ordered_ring(sizes5, order5, deg, src4, seg4$to, seg4$w)
          cost <- cost4 + ring_cost(blocks5, src4, seg4$to, seg4$w)
          if (cost < best_cost) {
            best_cost <- cost
            ring_orders <- list(`4`=order4, `5`=order5)
          }
        }
      }
      message("ring 4 order: ", paste(ring_orders[["4"]], collapse=", "))
      message("ring 5 order: ", paste(ring_orders[["5"]], collapse=", "))
    }
    rot <- seg %>% transmute(block=to, w_to, src=(a1 + a2) / 2)
    chain_rings[[k + 1]] <- ordered_ring(sizes, ring_orders[[as.character(k + 1)]], deg, rot$src, rot$block, rot$w_to)
    span_5 <- k + 1 == 5 && !is.null(ring5_span)
    if (span_5) {
      ## own scale: the blocks in the fixed order fill the span (half a gap at both ends)
      sizes5 <- sizes[match(ring_orders[["5"]], sizes$block), ]
      deg5 <- (diff(ring5_span) - nrow(sizes5) * block_gap) / sum(sizes5$size)
      chain_rings[[k + 1]] <- rotate_blocks(layout_blocks(setNames(sizes5$size, sizes5$block), deg5), ring5_span[1])
    }
    extra <- ring_extra_rotation[as.character(k + 1)]
    if (!is.na(extra) && !span_5) {
      chain_rings[[k + 1]] <- rotate_blocks(chain_rings[[k + 1]], extra)
    }
    chain_rings[[k + 1]] <- chain_rings[[k + 1]] %>%
      mutate(
        shift=coalesce(unname(block_extra_rotation[block]), 0),
        start=start + shift,
        end=end + shift
      ) %>%
      select(-shift)
  }
  chain_rings[[k + 1]] <- apply_block_start_fixed(chain_rings[[k + 1]])

  segments[[k]] <- seg %>%
    mutate(from_mid=(a1 + a2) / 2) %>%
    allocate_side(chain_rings[[k + 1]], "to", "from_mid", "b") %>%
    select(-from_center, -to_center_guess, -from_mid)
}

## current block positions, as start values for block_start_fixed (degrees clockwise from 12:00)
message("block positions (start - end, degrees clockwise from 12:00):")
iwalk(chain_rings[-1], function(blocks, i) {
  walk(seq_len(nrow(blocks)), function(j) {
    message(sprintf("  ring %d  %-22s %6.1f - %6.1f", i + 1, blocks$block[j], blocks$start[j] %% 360, blocks$end[j] %% 360))
  })
})

## flows: chains of segments through pass-through blocks, drawn as one continuous flow
all_segments <- bind_rows(segments)

chain_flow_polygons <- all_segments %>%
  filter(!is_pass(from) | seg == 1) %>%
  mutate(chain_id=row_number()) %>%
  pmap_dfr(function(...) {
    flow <- list(...)
    radii <- c(chain_radii[[flow$seg]][2], chain_radii[[flow$seg + 1]][1])
    left <- c(flow$a1, flow$b1)
    right <- c(flow$a2, flow$b2)
    current <- flow
    while (is_pass(current$to)) {
      nxt <- all_segments %>% filter(seg == current$seg + 1, from == current$to, site == current$site)
      if (nrow(nxt) == 0) break
      nxt <- as.list(nxt[1, ])
      radii <- c(radii, chain_radii[[nxt$seg]][2], chain_radii[[nxt$seg + 1]][1])
      left <- c(left, nxt$a1, nxt$b1)
      right <- c(right, nxt$a2, nxt$b2)
      current <- nxt
    }
    ## unwrap the angles along the chain (keeps the direction around the circle)
    for (i in seq_along(left)[-1]) {
      shift <- (left[i - 1] + right[i - 1]) / 2 + ang_diff((left[i - 1] + right[i - 1]) / 2, (left[i] + right[i]) / 2) -
        (left[i] + right[i]) / 2
      left[i] <- left[i] + shift
      right[i] <- right[i] + shift
    }
    ## one smooth curve from the first to the last ring: the pass-through slots only keep space free,
    ## the flow is not forced through them (that gave kinks)
    keep <- c(1, length(radii))
    radii <- radii[keep]
    left <- left[keep]
    right <- right[keep]
    edge <- function(angles) {
      t <- seq(0, 1, length.out=30)
      s <- 3 * t^2 - 2 * t^3
      map_dfr(seq_len(length(radii) - 1), function(i) {
        tibble(r=radii[i] + (radii[i + 1] - radii[i]) * t, angle=angles[i] + (angles[i + 1] - angles[i]) * s)
      })
    }
    l <- edge(left)
    r <- edge(right)
    n <- length(radii)
    bind_rows(
      to_xy(l$angle, l$r),
      arc_xy(left[n], right[n], radii[n]),
      to_xy(rev(r$angle), rev(r$r)),
      arc_xy(right[1], left[1], radii[1])
    ) %>%
      mutate(
        id=paste("chain", flow$chain_id),
        fill=case_when(
          flow$site %in% names(site_colors) ~ site_colors[flow$site],
          color_flows_by_source && flow$from %in% names(chain_block_colors) ~
            lighten(chain_block_colors[flow$from], outer_flow_lighten),
          TRUE ~ outer_flow_color
        )
      )
  })

chain_blocks <- imap_dfr(chain_rings, function(blocks, k) {
  blocks %>%
    filter(!is_pass(block)) %>%
    mutate(
      ring=k,
      fill=if (k == 1) site_colors[block] else chain_block_colors[block],
      r_inner=chain_radii[[k]][1],
      r_outer=chain_radii[[k]][2]
    )
}) %>%
  mutate(id=paste("block", row_number()))

chain_block_polygons <- chain_blocks %>%
  rowwise() %>%
  reframe(id, fill, block_poly(start, end, c(r_inner, r_outer)))

## labels: inside wide blocks (along the circle), next to narrow ones (radial, just outside)
is_dark <- function(colors) {
  rgb <- col2rgb(colors)
  (0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]) < 150
}

chain_labels <- chain_blocks %>%
  mutate(
    mid=((start + end) / 2) %% 360,
    wide=(end - start) >= if_else(ring == 1, 20, 14),
    unit=case_when(
      ring == 1 ~ " events",
      ring == 2 ~ " samples",
      ring == 3 ~ " blocks",
      ring == 5 ~ " cells",
      TRUE ~ " images"
    ),
    label=if_else(
      wide,
      paste0(block, "\n", prettyNum(round(size), big.mark=","), unit),
      paste0(block, " (", prettyNum(round(size), big.mark=","), ")")
    ),
    r_text=if_else(wide, (r_inner + r_outer) / 2, r_outer + 0.012),
    x=to_xy(mid, r_text)$x,
    y=to_xy(mid, r_text)$y,
    angle=if_else(wide, tangential_angle(mid), radial_angle(mid)),
    hjust=if_else(wide, 0.5, if_else(mid < 180, 0, 1)),
    color=if_else(wide & is_dark(fill), "white", "black")
  )

chain_legend <- c(if (chain_legend_sites) palette_legend(site_palette), palette_legend(chain_palette))
chain_legend <- chain_legend[!duplicated(names(chain_legend))]

chain_plot <- ggplot() +
  geom_polygon(data=chain_flow_polygons, aes(x, y, group=id, fill=fill), alpha=0.5, color=NA) +
  geom_polygon(data=chain_block_polygons, aes(x, y, group=id, fill=fill), color="white", linewidth=0.2) +
  geom_text(
    data=chain_labels,
    aes(x, y, label=label, angle=angle, hjust=hjust, color=color),
    size=1.8, lineheight=0.85
  ) +
  annotate(
    "text", 0, 0, size=2.2, lineheight=0.95,
    label=paste0(
      sum(ring1$size), " events\n",
      format(sum(chain_samples$size), big.mark=","), " samples\n",
      format(nrow(chain_images), big.mark=","), " images"
    )
  ) +
  ## legend from chain_palette (+ site_palette): the fill values are colours, the legend lists the
  ## palette colours with their labels (lightened flow colours are not listed)
  scale_fill_identity(
    guide=if (chain_legend_position == "none" || length(chain_legend) == 0) "none" else "legend",
    breaks=names(chain_legend),
    labels=unname(chain_legend),
    name=chain_legend_title
  ) +
  scale_color_identity() +
  ## room for the labels next to narrow outer blocks
  coord_equal(
    xlim=c(-1, 1) * (chain_radii[[5]][2] + 0.18),
    ylim=c(-1, 1) * (chain_radii[[5]][2] + 0.18)
  ) +
  theme_void() +
  theme(
    plot.background=element_rect(fill="white", color=NA),
    legend.position=chain_legend_position,
    legend.key.size=unit(3.5, "mm"),
    legend.text=element_text(size=7),
    legend.title=element_text(size=8)
  )

## a legend on the right / left needs extra width, so that the circle keeps its size
chain_width <- 7 + if (chain_legend_position %in% c("right", "left") && length(chain_legend) > 0) chain_legend_width else 0
ggsave(file.path(plot_dir, "biobank_overview_chain.png"), chain_plot, width=chain_width, height=7, dpi=400, bg="white")
ggsave_pdf(file.path(plot_dir, "biobank_overview_chain.pdf"), chain_plot, width=chain_width, height=7)


## ---- E: radial chain split by site: the rings of D, every ring split by site as in C -------------
## ring 1: sites (events), ring 2: sample types (samples), ring 3: resin blocks (samples),
## ring 4: imaging of the sample (images), ring 5: cells (images); same paths as D (segments_site).
## Every ring is split by site; the blocks of one site sit together, the site group is centred over its
## sources in the ring below (then pushed apart where groups overlap). Inside a site group the blocks are
## ordered by the position of their sources (ring 2: e_ring2_order). Paths that skip a ring run through
## an invisible pass-through gap of that ring. Rings 2 / 3 share one scale (degrees per sample),
## rings 4 / 5 one scale (degrees per image); the larger ring of each pair fills the circle.
## Flows are coloured by site.

## ---- settings E ----
e_site_gap <- 3       # gap between site groups (degrees)
e_node_gap <- 0.6     # gap between the blocks of one site (degrees)
## ring radii: site ring from e_inner_radius, all rings separated by the same gap, so that the outer
## edge of ring 5 is at e_outer_radius
e_inner_radius <- 0.15
e_outer_radius <- 0.71
e_ring_widths <- c(0.06, 0.05, 0.05, 0.05, 0.05)   # rings 1-5
e_ring_gap <- (e_outer_radius - e_inner_radius - sum(e_ring_widths)) / (length(e_ring_widths) - 1)
e_radii <- map(seq_along(e_ring_widths), function(k) {
  r0 <- e_inner_radius + sum(head(e_ring_widths, k - 1)) + (k - 1) * e_ring_gap
  c(r0, r0 + e_ring_widths[k])
})
## ring 3: one derived sample per MetaB sample (block name = key in radial_chain_palette)
## one derived sample per MetaB sample (block name = key in radial_chain_palette), in ring 3 or 4
## (ring 4: the flow passes ring 3 in an invisible gap)
e_metab_sequencing <- TRUE
e_metab_sequencing_block <- "V4V5 rRNA sequencing"
e_metab_sequencing_ring <- 4
## ring 2: plunge frozen split into cryo PF and cryo WM (FALSE: one PF block as in D)
e_split_pf <- TRUE
## order of the sample types inside a site (NULL: by source position)
e_ring2_order <- if (e_split_pf) c("MetaB", "HPF", "cryo PF", "cryo WM", "Chem fix") else ring2_order_fixed
e_flow_alpha <- 0.5
e_label_min <- c(`2`=4, `3`=4, `4`=4, `5`=4)  # rings 2-5: count inside blocks at least this wide (degrees)
e_site_label_min <- 20                 # ring 1: name + event count inside sites at least this wide, else inside the ring (towards the centre)
e_site_label_count <- TRUE             # ring 1: show the number of events under the site name
e_margin <- 0.04                       # space around the outer ring (no labels outside any more)
e_legend_position <- radial_chain_legend_position
## colours and legend: radial_chain_palette / radial_chain_ring_titles in the colour section at the top
e_block_colors <- palette_colors(radial_chain_palette)
radial_chain_palette <- radial_chain_palette %>%
  mutate(ring=if_else(key == e_metab_sequencing_block, e_metab_sequencing_ring, ring)) %>%
  filter(!(key == "PF" & e_split_pf), !(key %in% c("cryo PF", "cryo WM") & !e_split_pf))
## legend blocks: list of (title, entries with color + label), sites first if wanted, then ring 2-5
e_legend_blocks <- c(
  if (radial_chain_legend_sites) list(list(title=radial_chain_sites_title, pal=site_palette)),
  map(names(radial_chain_ring_titles), function(r) {
    list(title=radial_chain_ring_titles[[r]], pal=filter(radial_chain_palette, ring == as.numeric(r)))
  })
) %>%
  map(function(b) {
    b$pal <- b$pal %>% filter(!is.na(label)) %>% distinct(color, .keep_all=TRUE)
    if (!"repository" %in% names(b$pal)) b$pal$repository <- NA_character_
    b
  }) %>%
  keep(~nrow(.x$pal) > 0)

## ---- layout E ----
e_segments <- map(segments_site, ~mutate(filter(.x, site %in% site_levels), site=as.character(site)))

## plunge frozen split into cryo PF and cryo WM: samples per type (ring 2) and their cryo LM images
## (through the pass-through of ring 3 to cryo LM in ring 4), from the unmerged sample layer
if (e_split_pf) {
  pf_types <- c(`cryo PF`="cryo PF", `cryo WM`="cryo WM")
  pf_samples <- samples_per_type %>%
    filter(sample_type %in% pf_types, site %in% site_levels) %>%
    transmute(site=as.character(site), from=site, to=sample_type, w_from=as.numeric(size), w_to=as.numeric(size), seg=1L)
  pf_images <- imap(pf_types, ~images_of("cryo LM") %>% filter(sample_layer == .x, site %in% site_levels))
  e_segments[[1]] <- bind_rows(filter(e_segments[[1]], to != "PF"), pf_samples)
  e_segments[[2]] <- bind_rows(
    filter(e_segments[[2]], from != "PF"),
    imap_dfr(pf_images, ~per_site(.x, .y, paste(pass_prefix, .y), n_keys, n_keys))
  ) %>% mutate(seg=2L)
  e_segments[[3]] <- bind_rows(
    filter(e_segments[[3]], from != paste(pass_prefix, "PF")),
    imap_dfr(pf_images, ~per_site(.x, paste(pass_prefix, .y), "cryo LM", n_keys, n_img))
  ) %>% mutate(seg=3L)
  e_segments[1:3] <- map(e_segments[1:3], ~mutate(filter(.x, w_from > 0, w_to > 0), site=as.character(site), w_from=as.numeric(w_from), w_to=as.numeric(w_to)))
}

## derived samples that are not in D: one sequencing library per MetaB sample, in ring 3, or in ring 4
## through a pass-through of ring 3
if (e_metab_sequencing) {
  metab_n <- e_segments[[1]] %>% filter(to == "MetaB") %>% select(site, n=w_to)
  if (e_metab_sequencing_ring == 3) {
    e_segments[[2]] <- bind_rows(e_segments[[2]], transmute(metab_n, site, from="MetaB", to=e_metab_sequencing_block, w_from=n, w_to=n, seg=2L))
  } else {
    seq_pass <- paste(pass_prefix, "MetaB sequencing")
    e_segments[[2]] <- bind_rows(e_segments[[2]], transmute(metab_n, site, from="MetaB", to=seq_pass, w_from=n, w_to=n, seg=2L))
    e_segments[[3]] <- bind_rows(e_segments[[3]], transmute(metab_n, site, from=seq_pass, to=e_metab_sequencing_block, w_from=n, w_to=n, seg=3L))
  }
}

## ring sizes from the flows ending there
e_sizes <- map(e_segments, ~summarise(group_by(.x, site, block=to), size=sum(w_to), .groups="drop"))
e_gaps <- function(sizes) n_distinct(sizes$site) * e_site_gap + (nrow(sizes) - n_distinct(sizes$site)) * e_node_gap
e_deg <- function(a, b) min((360 - e_gaps(a)) / sum(a$size), (360 - e_gaps(b)) / sum(b$size))
e_deg_sample <- e_deg(e_sizes[[1]], e_sizes[[2]])
e_deg_image <- e_deg(e_sizes[[3]], e_sizes[[4]])

## centres pushed apart (circularly) until neighbours are at least (w_i + w_j) / 2 + gap apart
e_push_apart <- function(center, width, gap, iterations=5000) {
  o <- order(center %% 360)
  cen <- (center %% 360)[o]
  w <- width[o]
  n <- length(cen)
  if (n >= 2) {
    for (it in seq_len(iterations)) {
      moved <- FALSE
      for (i in seq_len(n)) {
        j <- if (i == n) 1 else i + 1
        need <- (w[i] + w[j]) / 2 + gap
        have <- (cen[j] - cen[i]) %% 360
        if (have < need - 1e-4) {
          push <- (need - have) / 2
          cen[i] <- cen[i] - push
          cen[j] <- cen[j] + push
          moved <- TRUE
        }
      }
      if (!moved) break
    }
  }
  out <- numeric(n)
  out[o] <- cen
  out
}

## blocks of a ring (site, block, size, center = centre of its sources) -> start / end
e_place_ring <- function(nodes, deg, fixed_order=NULL) {
  nodes <- nodes %>%
    mutate(width=size * deg) %>%
    group_by(site) %>%
    mutate(
      site_center=circular_center(center, size),
      rank=if (!is.null(fixed_order)) coalesce(match(block, fixed_order), 99L) + rank(ang_diff(site_center, center)) / 100
           else rank(ang_diff(site_center, center), ties.method="first")
    ) %>%
    arrange(site, rank) %>%
    ungroup()
  groups <- nodes %>%
    group_by(site) %>%
    summarise(center=first(site_center), width=sum(width) + (n() - 1) * e_node_gap, .groups="drop") %>%
    mutate(center=e_push_apart(center, width, e_site_gap), g_start=center - width / 2)
  nodes %>%
    left_join(groups %>% select(site, g_start), by="site") %>%
    group_by(site) %>%
    mutate(
      start=g_start + lag(cumsum(width + e_node_gap), default=0),
      end=start + width
    ) %>%
    ungroup() %>%
    select(site, block, size, start, end)
}

## pieces of the blocks of one side for their flows, ordered by the angle of `other` (relative to the
## block centre); fill = the flows fill the block (site ring: events vs samples)
e_allocate <- function(seg, ring, node_col, prefix, fill=FALSE) {
  width_col <- if (prefix == "a") "w_from" else "w_to"
  seg %>%
    left_join(ring %>% select(site, block, n_start=start, n_end=end, n_size=size), by=c("site", setNames("block", node_col))) %>%
    mutate(rel=ang_diff((n_start + n_end) / 2, other)) %>%
    arrange(site, .data[[node_col]], rel, to) %>%
    group_by(across(all_of(c("site", node_col)))) %>%
    group_modify(function(f, key) {
      size <- if (fill) sum(f[[width_col]]) else f$n_size[1]
      pieces <- split_node(f$n_start[1], f$n_end[1], size, f[[width_col]])
      f[[paste0(prefix, "1")]] <- pieces$from
      f[[paste0(prefix, "2")]] <- pieces$to
      f
    }) %>%
    ungroup() %>%
    select(-n_start, -n_end, -n_size, -rel)
}

e_rings <- list(ring1 %>% transmute(site=as.character(site), block=as.character(site), size, start, end))
for (k in 1:4) {
  seg <- e_segments[[k]]
  src <- e_rings[[k]]
  ## pass 1: source side in a provisional order -> centres of the target blocks -> target ring
  provisional <- seg %>%
    left_join(src %>% transmute(site, from=block, other=(start + end) / 2), by=c("site", "from")) %>%
    e_allocate(src, "from", "a", fill=(k == 1))
  nodes <- provisional %>%
    mutate(mid=(a1 + a2) / 2) %>%
    group_by(site, block=to) %>%
    summarise(size=sum(w_to), center=circular_center(mid, w_to), .groups="drop")
  e_rings[[k + 1]] <- e_place_ring(nodes, if (k <= 2) e_deg_sample else e_deg_image, if (k == 1) e_ring2_order)
  ## pass 2: source side ordered by the target blocks, target side by the source pieces
  e_segments[[k]] <- seg %>%
    left_join(e_rings[[k + 1]] %>% transmute(site, to=block, other=(start + end) / 2), by=c("site", "to")) %>%
    e_allocate(src, "from", "a", fill=(k == 1)) %>%
    mutate(other=(a1 + a2) / 2) %>%
    e_allocate(e_rings[[k + 1]], "to", "b") %>%
    select(-other)
}

## ---- drawing E ----
e_all <- bind_rows(e_segments)

## flows: from a visible block through pass-through gaps to the next visible block, one smooth curve
e_flow_polygons <- e_all %>%
  filter(!is_pass(from) | seg == 1) %>%
  mutate(flow_id=row_number()) %>%
  pmap_dfr(function(...) {
    flow <- list(...)
    r_from <- e_radii[[flow$seg]][2]
    left <- c(flow$a1, flow$b1)
    right <- c(flow$a2, flow$b2)
    current <- flow
    while (is_pass(current$to)) {
      nxt <- e_all %>% filter(seg == current$seg + 1, from == current$to, site == current$site)
      if (nrow(nxt) == 0) break
      nxt <- as.list(nxt[1, ])
      left <- c(left, nxt$a1, nxt$b1)
      right <- c(right, nxt$a2, nxt$b2)
      current <- nxt
    }
    r_to <- e_radii[[current$seg + 1]][1]
    ## first and last position only; the shorter way around the circle
    left <- left[c(1, length(left))]
    right <- right[c(1, length(right))]
    shift <- (left[1] + right[1]) / 2 + ang_diff((left[1] + right[1]) / 2, (left[2] + right[2]) / 2) - (left[2] + right[2]) / 2
    left[2] <- left[2] + shift
    right[2] <- right[2] + shift
    t <- seq(0, 1, length.out=40)
    s <- 3 * t^2 - 2 * t^3
    r <- r_from + (r_to - r_from) * t
    bind_rows(
      to_xy(left[1] + (left[2] - left[1]) * s, r),
      arc_xy(left[2], right[2], r_to),
      to_xy(rev(right[1] + (right[2] - right[1]) * s), rev(r)),
      arc_xy(right[1], left[1], r_from)
    ) %>%
      mutate(id=paste("e", flow$flow_id), fill=site_colors[flow$site])
  })

e_blocks <- imap_dfr(e_rings, function(blocks, k) {
  blocks %>%
    filter(!is_pass(block)) %>%
    mutate(
      ring=k,
      fill=if (k == 1) site_colors[block] else e_block_colors[block],
      r_inner=e_radii[[k]][1],
      r_outer=e_radii[[k]][2]
    )
}) %>%
  mutate(id=paste("eb", row_number()))

e_block_polygons <- e_blocks %>%
  rowwise() %>%
  reframe(id, fill, block_poly(start, end, c(r_inner, r_outer)))

e_labels <- bind_rows(
  ## sites: name (+ number of events) inside wide blocks (along the circle), "name (n)" inside the ring
  ## next to narrow ones
  e_blocks %>%
    filter(ring == 1) %>%
    mutate(
      wide=end - start >= e_site_label_min,
      label=case_when(
        !e_site_label_count ~ block,
        wide ~ paste0(block, "\n", size),
        TRUE ~ paste0(block, " (", size, ")")
      ),
      ## narrow sites: label on the inner side of the ring, pointing to the centre
      r_text=if_else(wide, (r_inner + r_outer) / 2, r_inner - 0.012), radial=!wide, inward=!wide
    ) %>%
    select(-wide),
  ## rings 2-5: count inside wide blocks
  e_blocks %>%
    filter(ring %in% 2:5, end - start >= e_label_min[as.character(ring)]) %>%
    mutate(label=prettyNum(round(size), big.mark=","), r_text=(r_inner + r_outer) / 2, radial=FALSE)
) %>%
  mutate(
    mid=((start + end) / 2) %% 360,
    x=to_xy(mid, r_text)$x,
    y=to_xy(mid, r_text)$y,
    angle=if_else(radial, radial_angle(mid), tangential_angle(mid)),
    inward=coalesce(inward, FALSE),
    hjust=case_when(!radial ~ 0.5, inward ~ if_else(mid < 180, 1, 0), TRUE ~ if_else(mid < 180, 0, 1)),
    color=if_else(!radial & is_dark(fill), "white", "black"),
    size=if_else(ring == 1, 1.8, 1.6),
    fontface="plain"
  )

e_plot <- ggplot() +
  geom_polygon(data=e_flow_polygons, aes(x, y, group=id, fill=fill), alpha=e_flow_alpha, color=NA) +
  geom_polygon(data=e_block_polygons, aes(x, y, group=id, fill=fill), color="white", linewidth=0.15) +
  geom_text(
    data=e_labels,
    aes(x, y, label=label, angle=angle, hjust=hjust, color=color, size=size, fontface=fontface),
    lineheight=0.85
  ) +
  annotate(
    "text", 0, 0, size=2.2, lineheight=0.95,
    label=paste0(
      sum(ring1$size), " events\n",
      format(sum(e_sizes[[1]]$size), big.mark=","), " samples\n",
      format(nrow(chain_images), big.mark=","), " images"
    )
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size_identity() +
  coord_equal(
    xlim=c(-1, 1) * (e_radii[[5]][2] + e_margin),
    ylim=c(-1, 1) * (e_radii[[5]][2] + e_margin)
  ) +
  theme_void() +
  theme(
    plot.background=element_rect(fill="white", color=NA),
    legend.position="none",
    legend.key.size=unit(3.5, "mm"),
    legend.text=element_text(size=7),
    legend.title=element_text(size=8, face="bold"),
    legend.box="vertical",
    legend.box.just="left",
    legend.justification="left",
    legend.spacing.y=unit(3, "mm")
  )

## ---- legend E: drawn as its own panel next to the circle ----
## per entry: colour box, label and the repository symbol(s) in an aligned column; below the blocks the
## repository key. Coordinates in inches (the panel is as high as the figure, 7 in).
e_fig_height <- 7
e_text_pt <- 7
e_title_pt <- 8
e_row_h <- 0.17          # height of one legend row (in)
e_block_gap <- 0.12      # space between the legend blocks (in)
e_symbol_gap <- 0.15     # distance between two repository symbols of one entry (in)

## width of texts in inches (measured on an off-screen device)
text_width_in <- function(x, size_pt, bold=FALSE) {
  if (length(x) == 0) return(0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  vapply(x, function(s) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(s, gp=grid::gpar(fontsize=size_pt, fontface=if (bold) "bold" else "plain"))),
      "inches", valueOnly=TRUE
    )
  }, numeric(1))
}

e_legend_rows <- imap_dfr(e_legend_blocks, function(b, i) {
  bind_rows(
    tibble(block=i, kind="title", label=b$title),
    b$pal %>% transmute(block=i, kind="entry", label, color, repository)
  )
})
e_repositories <- repository_palette %>%
  filter(repository %in% unlist(str_split(na.omit(e_legend_rows$repository), "\\s*;\\s*")))
e_legend_rows <- bind_rows(
  e_legend_rows,
  if (nrow(e_repositories) > 0) bind_rows(
    tibble(block=length(e_legend_blocks) + 1, kind="title", label=repository_legend_title),
    e_repositories %>% transmute(block=length(e_legend_blocks) + 1, kind="repository", label, repository, shape, color)
  )
) %>%
  mutate(
    row=row_number(),
    ## y from the top: one row per line, a gap before every block but the first
    y_top=(row - 1) * e_row_h + (block - 1) * e_block_gap
  )

x_box <- c(0, 0.13)
x_label <- 0.2
label_w <- max(text_width_in(e_legend_rows$label[e_legend_rows$kind == "entry"], e_text_pt), 0)
x_symbol <- x_label + label_w * 1.08 + 0.15   # measured widths are a bit smaller than in the png
max_symbols <- max(c(0, lengths(str_split(na.omit(e_legend_rows$repository[e_legend_rows$kind == "entry"]), ";"))))
e_legend_width <- max(
  x_symbol + max(max_symbols - 1, 0) * e_symbol_gap + 0.15,
  x_label + max(text_width_in(e_legend_rows$label[e_legend_rows$kind == "title"], e_title_pt, bold=TRUE), 0) + 0.05,
  x_label + max(text_width_in(e_legend_rows$label[e_legend_rows$kind == "repository"], e_text_pt), 0) + 0.1
)
## vertically centred
legend_h <- max(e_legend_rows$y_top) + e_row_h
e_legend_rows <- e_legend_rows %>% mutate(y=(e_fig_height + legend_h) / 2 - y_top - e_row_h / 2)

legend_symbols <- bind_rows(
  ## entries: symbols in the repository column
  e_legend_rows %>%
    filter(kind == "entry", !is.na(repository)) %>%
    select(y, repository_list=repository) %>%
    mutate(repository=str_split(repository_list, "\\s*;\\s*")) %>%
    unnest(repository) %>%
    group_by(y) %>%
    mutate(x=x_symbol + (row_number() - 1) * e_symbol_gap) %>%
    ungroup() %>%
    inner_join(repository_palette %>% select(repository, shape, color), by="repository"),
  ## repository key: symbol in the place of the colour box
  e_legend_rows %>% filter(kind == "repository") %>% transmute(y, repository, x=mean(x_box), shape, color)
)
unknown_repositories <- setdiff(
  unlist(str_split(na.omit(e_legend_rows$repository[e_legend_rows$kind == "entry"]), "\\s*;\\s*")),
  repository_palette$repository
)
if (length(unknown_repositories) > 0) warning("repositories not in repository_palette: ", paste(unknown_repositories, collapse=", "))

e_legend_plot <- ggplot() +
  geom_rect(
    data=filter(e_legend_rows, kind == "entry"),
    aes(xmin=x_box[1], xmax=x_box[2], ymin=y - e_row_h * 0.36, ymax=y + e_row_h * 0.36, fill=color)
  ) +
  geom_text(
    data=filter(e_legend_rows, kind == "title"),
    aes(x=x_box[1], y=y, label=label), hjust=0, size=e_title_pt / .pt, fontface="bold"
  ) +
  geom_text(
    data=filter(e_legend_rows, kind != "title"),
    aes(x=x_label, y=y, label=label), hjust=0, size=e_text_pt / .pt
  ) +
  geom_point(data=legend_symbols, aes(x, y, shape=shape, fill=color), color="white", size=2.4, stroke=0.3) +
  scale_fill_identity() +
  scale_shape_identity() +
  coord_cartesian(xlim=c(0, e_legend_width), ylim=c(0, e_fig_height), expand=FALSE) +
  theme_void() +
  theme(plot.margin=margin(0, 0, 0, 0), plot.background=element_rect(fill="white", color=NA))

if (e_legend_position != "none" && length(e_legend_blocks) > 0) {
  e_width <- e_fig_height + e_legend_width + 0.1
  e_figure <- patchwork::wrap_plots(e_plot, e_legend_plot, widths=c(e_fig_height, e_legend_width + 0.1))
} else {
  e_width <- e_fig_height
  e_figure <- e_plot
}

ggsave(file.path(plot_dir, "biobank_overview_radial_chain.png"), e_figure, width=e_width, height=e_fig_height, dpi=400, bg="white")
ggsave_pdf(file.path(plot_dir, "biobank_overview_radial_chain.pdf"), e_figure, width=e_width, height=e_fig_height)
