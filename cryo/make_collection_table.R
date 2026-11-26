library(tidyverse)
library(getopt)
library(jsonlite)

spec <- matrix(c(
  "all_s3", "d", 1, "character",
  "all_datasets", "a", 1, "character",
  "metadata_dir", "x", 1, "character",
  "s3_bucket", "b", 1, "character",
  "sheet_mode", "m", 1, "character",
  "google_key", "k", 1, "character",
  "collection_table_url", "u", 1, "character",
  "collection_table_sheet", "s", 1, "character",
  "local_collection_table", "l", 1, "character",
  "hpf_collection_table_sheet", "H", 1, "character",
  "local_hpf_collection_table", "L", 1, "character"
), ncol=4, byrow=TRUE)
opt <- getopt(spec)

value_or_default <- function(value, default) {
  if (is.null(value) || is.na(value) || value == "") default else value
}

sheet_mode <- value_or_default(opt$sheet_mode, "google")
collection_table_url <- value_or_default(opt$collection_table_url, "")
collection_table_sheet <- value_or_default(opt$collection_table_sheet, "cryo_collection_table")
local_collection_table <- value_or_default(opt$local_collection_table, "cryo_collection_table.tsv")
## HPF carrier screens (_HPF_ in the name) go into their own collection table
hpf_collection_table_sheet <- value_or_default(opt$hpf_collection_table_sheet, "cryo_hpf_collection_table")
local_hpf_collection_table <- value_or_default(opt$local_hpf_collection_table, "cryo_hpf_collection_table.tsv")
metadata_dir <- value_or_default(opt$metadata_dir, ".")
s3_bucket <- value_or_default(opt$s3_bucket, "s3test/imatrec/central_data_processing/cryo")

parse_mc_ls_path <- function(line) {
  parsed <- str_match(line, "^\\[.*?\\]\\s+\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2]
  ifelse(
    is.na(parsed),
    ## plain names (listing of the mounted S3 share) are used as they are
    coalesce(str_match(line, "^\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2], line),
    parsed
  )
}

## the mc alias (e.g. s3test) is dropped, the rest of the path is served by the public endpoint
s3_public_prefix <- function(bucket) {
  bucket_path <- bucket %>%
    str_remove("^[^/]+/") %>%
    str_remove("/$")
  file.path("https://schwab.buckets.embl.de", bucket_path)
}

metadata_row <- function(path) {
  metadata <- fromJSON(path, simplifyVector=FALSE)
  channels <- metadata$channels %||% list()
  channel_stats <- metadata$channel_stats %||% list()
  if (length(channels) > 0 && length(channel_stats) > 0) {
    stats_by_index <- set_names(channel_stats, map_chr(channel_stats, ~as.character(.x$index %||% "")))
    channels <- map(channels, function(channel) {
      stats <- stats_by_index[[as.character(channel$index %||% "")]]
      if (!is.null(stats)) {
        channel$min <- channel$min %||% stats$min
        channel$max <- channel$max %||% stats$max
        channel$contrast_limits <- channel$contrast_limits %||% stats$contrast_limits
      }
      channel
    })
  }
  tibble(
    name=metadata$name,
    metadata_raw_path=metadata$raw_path %||% "",
    x_scale_nm=metadata$x_scale_nm %||% NA_real_,
    y_scale_nm=metadata$y_scale_nm %||% NA_real_,
    z_scale_nm=metadata$z_scale_nm %||% NA_real_,
    shape=if (!is.null(metadata$shape)) paste(metadata$shape, collapse="x") else "",
    axes=metadata$axes %||% "",
    page_count=metadata$page_count %||% NA_integer_,
    source_suffix=metadata$source_suffix %||% "",
    size_c=metadata$size_c %||% NA_integer_,
    channels=list(channels)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

all_datasets <- read_tsv(opt$all_datasets, col_types=cols(.default=col_character())) %>%
  distinct(filename, .keep_all=TRUE)

s3_root_markers <- read_lines(opt$all_s3) %>%
  as_tibble() %>%
  mutate(
    s3_raw=parse_mc_ls_path(value),
    s3_raw=str_remove(s3_raw, "/$"),
    s3_prefix=str_match(s3_raw, "^([^/]+(?:\\.ome\\.zarr|\\.zarr))$")[, 2],
    name=s3_prefix %>%
      str_remove("\\.ome\\.zarr$") %>%
      str_remove("\\.zarr$")
  ) %>%
  filter(!is.na(name), name != "") %>%
  distinct(name, .keep_all=TRUE) %>%
  mutate(
    ## paste0, not file.path(..., "/"): on Linux that gives ".ome.zarr//", which the Fiji OME-Zarr reader cannot open
    uri=paste0(s3_public_prefix(s3_bucket), "/", s3_prefix, "/")
  ) %>%
  select(uri, name, s3_prefix, s3_raw)

## only datasets that are (still) in the input table; others on S3 (old test runs, removed rows) are left out
not_in_input <- setdiff(s3_root_markers$name, all_datasets$filename)
if (length(not_in_input) > 0) {
  message(
    "Leaving out ", length(not_in_input), " S3 dataset(s) that are not in the input table: ",
    paste(not_in_input, collapse=", ")
  )
}
s3_root_markers <- s3_root_markers %>%
  filter(name %in% all_datasets$filename)

if (nrow(s3_root_markers) == 0) {
  stop("No CRYO OME-Zarr datasets found in S3 listing that are in the input table.")
}

metadata_files <- list.files(
  metadata_dir,
  pattern="_metadata\\.json$",
  full.names=TRUE
)

metadata_table <- tibble(
  name=character(),
  metadata_raw_path=character(),
  x_scale_nm=double(),
  y_scale_nm=double(),
  z_scale_nm=double(),
  shape=character(),
  axes=character(),
  page_count=integer(),
  source_suffix=character(),
  size_c=integer(),
  channels=list()
)
if (length(metadata_files) > 0) {
  metadata_table <- metadata_files %>%
    map_dfr(metadata_row) %>%
    distinct(name, .keep_all=TRUE)
} else {
  warning("No *_metadata.json files found; collection table will use sheet metadata only.")
}

default_channel_colors <- c("red", "green", "yellow", "blue", "magenta", "cyan", "white")

near_square_columns <- function(n) {
  if (n <= 1) {
    return(1L)
  }
  as.integer(ceiling(sqrt(n)))
}

color_for_display <- function(display, fallback) {
  compact <- str_replace_all(tolower(display %||% ""), "[^a-z0-9]+", "")
  case_when(
    str_detect(compact, "gfp") ~ "cyan",
    str_detect(compact, "tl") ~ "white",
    str_detect(compact, "reflection") ~ "white",
    str_detect(compact, "chloa") | str_detect(compact, "chlorophyll") ~ "magenta",
    str_detect(compact, "568|594|rfp|mcherry") ~ "magenta",
    str_detect(compact, "pe") ~ "yellow",
    TRUE ~ fallback
  )
}

rounded_wavelength <- function(value) {
  value <- suppressWarnings(as.numeric(value %||% NA_real_))
  ifelse(is.na(value), NA_real_, round(value))
}

canonical_channel_display <- function(channel, fallback_display) {
  text <- paste(
    channel$display %||% "",
    channel$label %||% "",
    channel$fluor %||% "",
    sep=" "
  )
  compact <- str_replace_all(tolower(text), "[^a-z0-9]+", "")
  excitation <- rounded_wavelength(channel$excitation_wavelength_nm %||% NA_real_)
  emission <- rounded_wavelength(channel$emission_wavelength_nm %||% NA_real_)

  case_when(
    ## reflection channel of the HPF screens (Reflection-T2)
    str_detect(compact, "reflection") ~ "Reflection",
    str_detect(compact, "tl") & !str_detect(compact, "tpmt") ~ "TL",
    str_detect(compact, "gfp") | str_detect(compact, "egfp") |
      excitation == 488 | emission == 509 ~ "GFP",
    str_detect(compact, "alexafluor555") | str_detect(compact, "af555") |
      compact %in% c("pe") | excitation == 553 | emission == 568 ~ "PE",
    str_detect(compact, "chloa") | str_detect(compact, "chlorophylla") |
      excitation == 655 | emission == 667 ~ "ChloA",
    TRUE ~ fallback_display
  )
}

sanitize_channel_display <- function(value, index) {
  value <- value %||% paste0("channel_", index)
  value <- str_replace_all(as.character(value), "[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  ifelse(value == "", paste0("channel_", index), value)
}

normalize_channels <- function(channels, size_c=NA_integer_) {
  if (is.null(channels) || length(channels) == 0) {
    return(list(list(index=0L, label="channel_0", display="channel_0", color=default_channel_colors[[1]])))
  }

  normalized <- map(seq_along(channels), function(i) {
    channel <- channels[[i]]
    index <- as.integer(channel$index %||% (i - 1L))
    label <- channel$label %||% paste0("channel_", index)
    display <- channel$display %||% sanitize_channel_display(label, index)
    color <- channel$color %||% default_channel_colors[[(index %% length(default_channel_colors)) + 1]]
    contrast_limits <- channel$contrast_limits %||% ""
    display <- canonical_channel_display(channel, sanitize_channel_display(display, index))
    list(
      index=index,
      label=as.character(label),
      display=display,
      color=as.character(color_for_display(display, color)),
      contrast_limits=as.character(contrast_limits)
    )
  })

  ## `index` is the channel's position on the image C axis and is what MoBIE loads, so channels
  ## are neither dropped nor renumbered here (that used to shift channels onto the wrong data).
  normalized <- normalized[order(map_int(normalized, "index"))]
  if (!is.na(size_c) && size_c > 0) {
    normalized <- keep(normalized, ~.x$index < size_c)
  }
  ## two channels with the same display name (e.g. two TL tracks) get the index appended
  displays <- map_chr(normalized, "display")
  duplicated_display <- displays %in% displays[duplicated(tolower(displays))]
  normalized <- map2(normalized, duplicated_display, function(channel, is_duplicate) {
    if (is_duplicate) channel$display <- paste0(channel$display, "_", channel$index)
    channel
  })
  normalized
}

## optional input-table columns: site (MoBIE view) and grid (MoBIE grid, defaults to the site).
## Without them the site is taken from the leading letters of the name (e.g. Roscoff_HPF_...).
for (optional_column in c("site", "grid")) {
  if (!optional_column %in% names(all_datasets)) {
    all_datasets[[optional_column]] <- NA_character_
  }
}
all_datasets <- all_datasets %>%
  mutate(across(c(site, grid), ~na_if(str_trim(.x), "")))

base_table <- s3_root_markers %>%
  left_join(
    all_datasets %>%
      select(
        name=filename,
        shortname,
        raw_path,
        output_path,
        req_mem,
        sheet_x_scale=x_scale,
        sheet_y_scale=y_scale,
        sheet_z_scale=z_scale,
        sheet_site=site,
        sheet_grid=grid
      ),
    by="name"
  ) %>%
  left_join(metadata_table, by="name") %>%
  mutate(
    x_scale_nm=coalesce(as.character(x_scale_nm), sheet_x_scale),
    y_scale_nm=coalesce(as.character(y_scale_nm), sheet_y_scale),
    z_scale_nm=coalesce(as.character(z_scale_nm), sheet_z_scale),
    source_name=name,
    site=coalesce(sheet_site, str_extract(source_name, "^[A-Za-z]+")),
    collection=if_else(str_detect(source_name, "_HPF_"), "hpf", "pf_wm"),
    view=coalesce(site, "cryo"),
    grid=coalesce(sheet_grid, site, "cryo"),
    channels=map2(channels, size_c, normalize_channels),
    exclusive=FALSE,
    blend="sum",
    format="OmeZarr",
    type="intensities"
  ) %>%
  group_by(collection, grid) %>%
  mutate(
    ## position within the site's grid (per collection table)
    grid_index=row_number() - 1L,
    grid_columns=near_square_columns(n()),
    grid_position=paste0("(", grid_index %% grid_columns, ",", grid_index %/% grid_columns, ")")
  ) %>%
  ungroup() %>%
  select(-grid_columns)

col_table <- base_table %>%
  unnest_longer(channels) %>%
  mutate(
    channel=map_int(channels, "index"),
    channel_label=map_chr(channels, "label"),
    display=map_chr(channels, "display"),
    color=map_chr(channels, "color"),
    contrast_limits=map_chr(channels, "contrast_limits"),
    name=paste0(source_name, "_c", channel, "_", display),
    exclusive=T
  ) %>%
  select(-channels) %>%
  select(
    collection,
    uri,
    name,
    view,
    grid,
    grid_position,
    channel,
    display,
    color,
    contrast_limits,
    blend,
    format,
    exclusive
  )

pf_wm_table <- col_table %>% filter(collection == "pf_wm") %>% select(-collection)
hpf_table <- col_table %>% filter(collection == "hpf") %>% select(-collection)
message("Collection tables: ", n_distinct(pf_wm_table$uri), " PF/WM and ", n_distinct(hpf_table$uri), " HPF dataset(s)")

write_tsv(pf_wm_table, file=local_collection_table)
write_tsv(hpf_table, file=local_hpf_collection_table)

if (sheet_mode == "google") {
  if (collection_table_url == "") {
    stop("--collection_table_url is required when --sheet_mode google")
  }
  json_key <- opt$google_key
  if (is.null(json_key) || is.na(json_key)) {
    stop("--google_key is required when --sheet_mode google")
  }

  library(googlesheets4)
  library(googledrive)
  gs4_auth(path=json_key)
  drive_auth(path=json_key)
  ## an empty table would replace the tab with a header-only sheet, so only non-empty ones are written
  if (nrow(pf_wm_table) > 0) {
    write_sheet(pf_wm_table, ss=collection_table_url, sheet=collection_table_sheet)
  }
  if (nrow(hpf_table) > 0) {
    write_sheet(hpf_table, ss=collection_table_url, sheet=hpf_collection_table_sheet)
  }
}

write_tsv(tibble(done="done"), file="done.tsv")
