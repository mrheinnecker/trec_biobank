

library(tidyverse)
library(circlize)
library(ggalluvial)

## inputs written by biobank_overview.R
ov_dir <- "C:/projects/tem_classification/ov"
plot_dir <- "C:/projects/tem_classification/plots"
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

## a PDF that is open in a viewer cannot be overwritten on Windows: warn instead of stopping
ggsave_pdf <- function(filename, plot, ...) {
  tryCatch(
    ggsave(filename, plot, ...),
    error=function(e) warning(basename(filename), " not written (open in a viewer?)")
  )
}

## site order and colors as in cryo/cryo_ov_figure.R
site_colors <- c(
  Roscoff="#332288",
  Tallinn="#4477AA",
  Kristineberg="#66CCEE",
  Bilbao="#228833",
  Porto="#AA9933",
  Barcelona="#EE7733",
  Naples="#CC6677",
  Athens="#AA3377"
)

coverage <- read_tsv(file.path(ov_dir, "event_coverage.tsv"), show_col_types=FALSE)
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
sample_col_fun <- colorRamp2(count_breaks, c("grey94", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))
image_col_fun <- colorRamp2(count_breaks, c("grey94", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))

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

sample_type_colors <- c(
  HPF="#2166AC",
  `cryo PF`="#4393C3",
  `cryo WM`="#92C5DE",
  MetaB="#5AAE61",
  `Chem fix`="#9970AB"
)
imaging_colors <- c(
  TEM="#7F2704",
  HiTT="#D94801",
  `resin LM`="#FD8D3C",
  `cryo LM PF`="#FDAE6B",
  `cryo LM WM`="#FDD0A2",
  `cryo LM HPF`="#FEE6CE",
  SEM="#8C6D31",
  ExM="#35978F"
)
event_only <- "event only"
event_only_color <- "grey75"

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
    breaks=c(sample_type_colors, event_only_color, imaging_colors),
    labels=c(names(sample_type_colors), "no sample record", names(imaging_colors)),
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

chain_block_colors <- c(
  HPF="#2166AC",
  PF="#4393C3",
  MetaB="#5AAE61",
  `Chem fix`="#9970AB",
  Epon="#B2182B",
  Lowicryl="#EF8A62",
  HiTT="#D94801",
  `resin LM`="#FD8D3C",
  `cryo LM`="#FDAE6B",
  TEM="#7F2704",
  SEM="#8C6D31",
  ExM="#35978F"
)
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
key_resin <- bind_rows(
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
    distinct(site, sample_key)
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
  scale_fill_identity() +
  scale_color_identity() +
  ## room for the labels next to narrow outer blocks
  coord_equal(
    xlim=c(-1, 1) * (chain_radii[[5]][2] + 0.18),
    ylim=c(-1, 1) * (chain_radii[[5]][2] + 0.18)
  ) +
  theme_void() +
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(plot_dir, "biobank_overview_chain.png"), chain_plot, width=7, height=7, dpi=400, bg="white")
ggsave_pdf(file.path(plot_dir, "biobank_overview_chain.pdf"), chain_plot, width=7, height=7)
