#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(googlesheets4)
  library(grid)
})

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      out$help <- TRUE
      i <- i + 1L
    } else if (arg %in% c("--anonymous", "-a")) {
      out$anonymous <- TRUE
      i <- i + 1L
    } else if (str_detect(arg, "^--[^=]+=")) {
      key <- str_remove(str_extract(arg, "^--[^=]+"), "^--") %>%
        str_replace_all("-", "_")
      out[[key]] <- str_remove(arg, "^--[^=]+=")
      i <- i + 1L
    } else if (str_detect(arg, "^--")) {
      key <- str_remove(arg, "^--") %>%
        str_replace_all("-", "_")
      if (i == length(args) || str_detect(args[[i + 1L]], "^--")) {
        stop("Missing value for ", arg)
      }
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      stop("Unknown positional argument: ", arg)
    }
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly=TRUE))

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript annotations/collection_table_overview.R \\\n",
    "    --google_key /path/to/service-account.json \\\n",
    "    --outdir annotations/collection_overview\n\n",
    "Options:\n",
    "  --google_key PATH   Google service-account JSON key. Defaults to GOOGLE_KEY,\n",
    "                      then ./trec-tem-screen-e98a2e03f58b.json if present.\n",
    "  --outdir PATH       Output directory. Default: annotations/collection_overview\n",
    "  --prefix VALUE      Output filename prefix. Default: collection_table_overview\n",
    "  --modalities CSV    Optional subset, e.g. TEM,HITT,CRYO.\n",
    "  --anonymous         Use unauthenticated googlesheets4 access.\n",
    "  --help              Show this message.\n",
    sep=""
  )
}

if (!is.null(opt$help) && isTRUE(opt$help)) {
  usage()
  quit(save="no", status=0)
}

value_or_default <- function(value, default) {
  if (is.null(value) || is.na(value) || identical(value, "")) default else value
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

script_arg <- commandArgs(trailingOnly=FALSE) %>%
  keep(~str_detect(.x, "^--file=")) %>%
  first()
script_path <- if (is.null(script_arg)) {
  file.path(getwd(), "annotations", "collection_table_overview.R")
} else {
  str_remove(script_arg, "^--file=")
}
script_path <- normalizePath(script_path, mustWork=FALSE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork=FALSE)
default_key <- file.path(repo_root, "trec-tem-screen-e98a2e03f58b.json")

google_key <- value_or_default(opt$google_key, Sys.getenv("GOOGLE_KEY", unset=""))
if (google_key == "" && file.exists(default_key)) {
  google_key <- default_key
}

outdir <- value_or_default(opt$outdir, file.path(repo_root, "annotations", "collection_overview"))
prefix <- value_or_default(opt$prefix, "collection_table_overview")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

collection_sources <- tribble(
  ~modality, ~sheet_url, ~sheet_name,
  "TEM", "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1643802951#gid=1643802951", "tem_collection_table",
  "SEM", "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1595862976#gid=1595862976", "sem_collection_table",
  "HITT", "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1062275333#gid=1062275333", "hitt_collection_table",
  "CRYO", "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=199938698#gid=199938698", "cryo_collection_table",
  "PLASTIC", "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=199938698#gid=199938698", "plastic_collection_table",
  "EXM", "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit", "exm_collection_table"
)

if (!is.null(opt$modalities) && !is.na(opt$modalities) && opt$modalities != "") {
  keep_modalities <- str_split(opt$modalities, ",")[[1]] %>%
    str_trim() %>%
    str_to_upper()
  collection_sources <- collection_sources %>%
    filter(modality %in% keep_modalities)
}

if (nrow(collection_sources) == 0) {
  stop("No collection sources selected.")
}

if (!is.null(opt$anonymous) && isTRUE(opt$anonymous)) {
  gs4_deauth()
} else {
  if (google_key == "" || !file.exists(google_key)) {
    stop(
      "No Google key found. Pass --google_key /path/to/key.json, set GOOGLE_KEY, ",
      "or run with --anonymous for public sheets."
    )
  }
  gs4_auth(path=google_key)
}

non_empty_rows <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }
  df %>%
    mutate(across(everything(), ~replace_na(as.character(.x), ""))) %>%
    filter(if_any(everything(), ~str_trim(.x) != ""))
}

first_existing_column <- function(df, candidates) {
  matches <- candidates[candidates %in% names(df)]
  if (length(matches) == 0) {
    NULL
  } else {
    matches[[1]]
  }
}

dataset_id_from_uri <- function(uri) {
  uri %>%
    as.character() %>%
    str_remove("[?#].*$") %>%
    str_remove("/+$") %>%
    basename() %>%
    str_remove("\\.ome\\.zarr$") %>%
    str_remove("\\.zarr$") %>%
    na_if("")
}

derive_dataset_id <- function(df) {
  source_col <- first_existing_column(
    df,
    c("source_name", "dataset", "dataset_name", "image_name", "filename", "source")
  )

  if (!is.null(source_col)) {
    dataset_id <- df[[source_col]] %>%
      as.character() %>%
      str_trim() %>%
      na_if("")
  } else if ("uri" %in% names(df)) {
    dataset_id <- dataset_id_from_uri(df$uri)
  } else if ("name" %in% names(df)) {
    dataset_id <- df$name %>%
      as.character() %>%
      str_trim() %>%
      str_remove("_c[0-9]+_[A-Za-z0-9]+$") %>%
      na_if("")
  } else {
    dataset_id <- paste0("row_", seq_len(nrow(df)))
  }

  if ("uri" %in% names(df)) {
    uri_id <- dataset_id_from_uri(df$uri)
    dataset_id <- coalesce(dataset_id, uri_id)
  }

  missing <- is.na(dataset_id) | dataset_id == ""
  dataset_id[missing] <- paste0("row_", which(missing))
  dataset_id
}

derive_site <- function(df, dataset_id) {
  if ("site" %in% names(df)) {
    site <- df$site %>%
      as.character() %>%
      str_trim() %>%
      na_if("")
  } else {
    site <- rep(NA_character_, length(dataset_id))
  }

  coalesce(
    site,
    str_extract(dataset_id, "ATH|BAR|BIL|KRI|NAP|POR|ROS|TAL|TES|VIG|Vigo") %>%
      str_to_upper() %>%
      str_replace("^VIGO$", "VIG")
  ) %>%
    replace_na("unknown")
}

read_collection_table <- function(modality, sheet_url, sheet_name) {
  message("Reading ", modality, " / ", sheet_name)
  raw <- read_sheet(sheet_url, sheet=sheet_name, col_types="c")
  df <- non_empty_rows(raw)
  dataset_id <- derive_dataset_id(df)
  site <- derive_site(df, dataset_id)

  df %>%
    mutate(
      modality=modality,
      dataset_id=dataset_id,
      site=site,
      sheet_name=sheet_name
    )
}

tables <- pmap(
  collection_sources,
  function(modality, sheet_url, sheet_name) {
    tryCatch(
      list(
        status="ok",
        modality=modality,
        table=read_collection_table(modality, sheet_url, sheet_name),
        error=NA_character_
      ),
      error=function(e) {
        warning("Failed to read ", modality, ": ", conditionMessage(e))
        list(
          status="failed",
          modality=modality,
          table=tibble(),
          error=conditionMessage(e)
        )
      }
    )
  }
)

read_status <- map_dfr(tables, ~tibble(
  modality=.x$modality,
  status=.x$status,
  error=.x$error
))

collection_rows <- map_dfr(tables, "table")

if (nrow(collection_rows) == 0) {
  write_tsv(read_status, file.path(outdir, paste0(prefix, "_read_status.tsv")))
  stop("No collection table rows could be read. See read-status TSV for errors.")
}

modalities_order <- collection_sources$modality
palette <- c(
  TEM="#3B6EA8",
  SEM="#7A4EAB",
  HITT="#C85A3D",
  CRYO="#2F8F83",
  PLASTIC="#D9A441",
  EXM="#4F7F3A"
)

pick_first_column <- function(df, columns) {
  values <- rep(NA_character_, nrow(df))
  for (column in columns) {
    if (column %in% names(df)) {
      column_values <- df[[column]] %>%
        as.character() %>%
        str_trim() %>%
        na_if("")
      values <- coalesce(values, column_values)
    }
  }
  values
}

clean_size_fraction <- function(value) {
  value <- value %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("-", "to") %>%
    str_replace_all("_", "to") %>%
    str_replace_all("(?i)um", "") %>%
    na_if("")
  coalesce(value, "unknown")
}

clean_sampling_time <- function(value) {
  value <- value %>%
    as.character() %>%
    str_trim() %>%
    str_to_upper() %>%
    str_replace_all("_", "") %>%
    na_if("")
  case_when(
    value %in% c("AM", "MORNING") ~ "AM",
    value %in% c("PM", "AFTERNOON") ~ "PM",
    value %in% c("MID", "MIDDAY", "NOON") ~ "MID",
    value %in% c("TARA") ~ "TARA",
    TRUE ~ coalesce(value, "unknown")
  )
}

collection_rows$size_fraction <- clean_size_fraction(
  pick_first_column(collection_rows, c("size_frac", "size_fraction", "Size fraction"))
)
collection_rows$sampling_time_group <- clean_sampling_time(
  pick_first_column(collection_rows, c("sampling_time", "Sampling time", "sampling"))
)

dataset_metadata <- collection_rows %>%
  distinct(modality, dataset_id, .keep_all=TRUE)

summary_counts <- collection_rows %>%
  group_by(modality) %>%
  summarise(
    collection_rows=n(),
    datasets=n_distinct(dataset_id),
    sites=n_distinct(site[site != "unknown"]),
    rows_per_dataset=collection_rows / datasets,
    .groups="drop"
  ) %>%
  mutate(
    modality=factor(modality, levels=modalities_order),
    label=scales::comma(datasets)
  ) %>%
  arrange(modality)

site_levels <- dataset_metadata %>%
  filter(site != "unknown") %>%
  distinct(site) %>%
  pull(site) %>%
  sort()

site_counts <- dataset_metadata %>%
  filter(site != "unknown") %>%
  count(modality, site, name="datasets") %>%
  complete(
    modality=modalities_order,
    site=site_levels,
    fill=list(datasets=0)
  ) %>%
  mutate(
    modality=factor(modality, levels=modalities_order),
    site=factor(site, levels=rev(site_levels))
  )

size_fraction_counts <- dataset_metadata %>%
  count(modality, size_fraction, name="datasets") %>%
  group_by(modality) %>%
  mutate(size_fraction_rank=rank(-datasets, ties.method="first")) %>%
  ungroup() %>%
  mutate(
    size_fraction=if_else(size_fraction_rank <= 8, size_fraction, "other"),
    modality=factor(modality, levels=modalities_order)
  ) %>%
  group_by(modality, size_fraction) %>%
  summarise(datasets=sum(datasets), .groups="drop")

sampling_time_counts <- dataset_metadata %>%
  count(modality, sampling_time_group, name="datasets") %>%
  group_by(modality) %>%
  mutate(sampling_time_rank=rank(-datasets, ties.method="first")) %>%
  ungroup() %>%
  mutate(
    sampling_time_group=if_else(sampling_time_rank <= 8, sampling_time_group, "other"),
    modality=factor(modality, levels=modalities_order)
  ) %>%
  group_by(modality, sampling_time_group) %>%
  summarise(datasets=sum(datasets), .groups="drop")

write_tsv(collection_rows, file.path(outdir, paste0(prefix, "_collection_rows.tsv")))
write_tsv(dataset_metadata, file.path(outdir, paste0(prefix, "_dataset_metadata.tsv")))
write_tsv(summary_counts, file.path(outdir, paste0(prefix, "_summary_counts.tsv")))
write_tsv(site_counts, file.path(outdir, paste0(prefix, "_site_counts.tsv")))
write_tsv(size_fraction_counts, file.path(outdir, paste0(prefix, "_size_fraction_counts.tsv")))
write_tsv(sampling_time_counts, file.path(outdir, paste0(prefix, "_sampling_time_counts.tsv")))
write_tsv(read_status, file.path(outdir, paste0(prefix, "_read_status.tsv")))

base_theme <- theme_minimal(base_size=13) +
  theme(
    panel.grid.minor=element_blank(),
    plot.title.position="plot",
    plot.title=element_text(face="bold", size=15),
    plot.subtitle=element_text(color="grey30"),
    axis.title=element_text(face="bold"),
    strip.text=element_text(face="bold"),
    legend.position="bottom"
  )

p_counts <- summary_counts %>%
  ggplot(aes(x=reorder(modality, datasets), y=datasets, fill=modality)) +
  geom_col(width=0.72, color="white", linewidth=0.5) +
  geom_text(aes(label=scales::comma(datasets)), hjust=-0.12, fontface="bold") +
  coord_flip(clip="off") +
  scale_fill_manual(values=palette, guide="none") +
  scale_y_continuous(labels=scales::comma, expand=expansion(mult=c(0, 0.18))) +
  labs(
    title="Processed datasets by modality",
    subtitle="One dataset is counted once, even when a collection table has one row per channel.",
    x=NULL,
    y="Unique datasets"
  ) +
  base_theme

p_sites <- site_counts %>%
  ggplot(aes(x=site, y=datasets, fill=modality)) +
  geom_col(width=0.72, color="white", linewidth=0.3) +
  geom_text(
    aes(label=if_else(datasets > 0, scales::comma(datasets), "")),
    hjust=-0.12,
    size=3.1,
    fontface="bold"
  ) +
  coord_flip(clip="off") +
  facet_wrap(vars(modality)) +
  scale_fill_manual(values=palette, guide="none") +
  scale_y_continuous(labels=scales::comma, expand=expansion(mult=c(0, 0.18))) +
  labs(
    title="Dataset composition by site",
    subtitle="Sites are read from the collection table when present, otherwise inferred from dataset names.",
    x=NULL,
    y="Unique datasets"
  ) +
  base_theme

p_size_fraction <- size_fraction_counts %>%
  mutate(size_fraction=fct_reorder(size_fraction, datasets, .fun=sum)) %>%
  ggplot(aes(x=size_fraction, y=datasets, fill=modality)) +
  geom_col(width=0.72, color="white", linewidth=0.3) +
  coord_flip() +
  facet_wrap(vars(modality), scales="free_y") +
  scale_fill_manual(values=palette, guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(
    title="Dataset composition by size fraction",
    subtitle="Uses size_frac / Size fraction columns when present; otherwise marked unknown.",
    x=NULL,
    y="Unique datasets"
  ) +
  base_theme

p_sampling_time <- sampling_time_counts %>%
  mutate(sampling_time_group=fct_relevel(sampling_time_group, "AM", "MID", "PM", "TARA", "unknown", after=0)) %>%
  ggplot(aes(x=sampling_time_group, y=datasets, fill=modality)) +
  geom_col(width=0.72, color="white", linewidth=0.3) +
  facet_wrap(vars(modality), scales="free_y") +
  scale_fill_manual(values=palette, guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(
    title="Dataset composition by sampling time",
    subtitle="Uses sampling_time-like collection-table columns when present; otherwise marked unknown.",
    x=NULL,
    y="Unique datasets"
  ) +
  base_theme +
  theme(axis.text.x=element_text(angle=25, hjust=1))

save_plot <- function(plot, suffix, width, height) {
  ggsave(
    filename=file.path(outdir, paste0(prefix, "_", suffix, ".png")),
    plot=plot,
    width=width,
    height=height,
    dpi=220
  )
  ggsave(
    filename=file.path(outdir, paste0(prefix, "_", suffix, ".pdf")),
    plot=plot,
    width=width,
    height=height,
    device=cairo_pdf
  )
}

save_plot(p_counts, "datasets_by_modality", 8.5, 5.5)
save_plot(p_sites, "datasets_by_site", 10.5, 7.0)
save_plot(p_size_fraction, "datasets_by_size_fraction", 10.5, 7.0)
save_plot(p_sampling_time, "datasets_by_sampling_time", 10.5, 6.5)

unlink(file.path(
  outdir,
  paste0(
    prefix,
    c(
      "_rows_vs_datasets.png", "_rows_vs_datasets.pdf",
      "_rows_per_dataset.png", "_rows_per_dataset.pdf"
    )
  )
))

combined_png <- file.path(outdir, paste0(prefix, "_overview.png"))
combined_pdf <- file.path(outdir, paste0(prefix, "_overview.pdf"))

draw_combined <- function(device_fun, filename) {
  device_fun(filename)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(
    nrow=1,
    ncol=2,
    widths=unit(c(0.95, 1.05), "null")
  )))
  print(p_counts, vp=viewport(layout.pos.row=1, layout.pos.col=1))
  print(p_sites, vp=viewport(layout.pos.row=1, layout.pos.col=2))
  dev.off()
}

draw_combined(
  function(filename) png(filename, width=3200, height=1300, res=220),
  combined_png
)
draw_combined(
  function(filename) cairo_pdf(filename, width=14.5, height=5.9),
  combined_pdf
)

message("Wrote:")
message("  ", file.path(outdir, paste0(prefix, "_summary_counts.tsv")))
message("  ", combined_png)
message("  ", combined_pdf)
