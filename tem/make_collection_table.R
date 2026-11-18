library(tidyverse)
library(getopt)

spec <- matrix(c(
  # long option                  short  arg  type
  "all_s3", "d",   1,   "character",
  "sheet_mode", "m", 1, "character",
  "google_key", "k", 1, "character",
  "collection_table_url", "u", 1, "character",
  "collection_table_sheet", "s", 1, "character",
  "local_collection_table", "l", 1, "character",
  "image_log_url", "i", 1, "character",
  "image_log_sheet", "g", 1, "character",
  "local_image_log", "a", 1, "character",
  "image_stats_dir", "x", 1, "character",
  "expected_datasets", "e", 1, "character",
  "s3_bucket", "b", 1, "character"
),
ncol = 4,
byrow = TRUE)
opt <- getopt(spec)

#opt$all_s3 <- "/g/schwab/marco/repos/trec_biobank/all_s3_entries.txt"
#opt$google_key <- "/g/schwab/marco/repos/trec_biobank/tem/trec-tem-screen-e98a2e03f58b.json"

sheet_mode <- opt$sheet_mode
if (is.null(sheet_mode) || is.na(sheet_mode)) {
  sheet_mode <- "local"
}

collection_table <- opt$collection_table_url
if (is.null(collection_table) || is.na(collection_table)) {
  collection_table <- "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1643802951#gid=1643802951"
}

collection_table_sheet <- opt$collection_table_sheet
if (is.null(collection_table_sheet) || is.na(collection_table_sheet)) {
  collection_table_sheet <- "collection_table"
}

local_collection_table <- opt$local_collection_table
if (is.null(local_collection_table) || is.na(local_collection_table)) {
  local_collection_table <- "collection_table.tsv"
}

image_log_url <- opt$image_log_url
if (is.null(image_log_url) || is.na(image_log_url)) {
  image_log_url <- "https://docs.google.com/spreadsheets/d/143uVeeJ72SQE5eK01lzWYCEiT7pJUF3lX7hJl3R9s9I/edit?gid=258669282#gid=258669282"
}

image_log_sheet <- opt$image_log_sheet
if (is.null(image_log_sheet) || is.na(image_log_sheet)) {
  image_log_sheet <- "image_log"
}

local_image_log <- opt$local_image_log
if (is.null(local_image_log) || is.na(local_image_log)) {
  local_image_log <- "image_log_local.tsv"
}

image_stats_dir <- opt$image_stats_dir
if (is.null(image_stats_dir) || is.na(image_stats_dir)) {
  image_stats_dir <- "."
}

s3_bucket <- opt$s3_bucket
if (is.null(s3_bucket) || is.na(s3_bucket)) {
  s3_bucket <- "s3test/imatrec/central_data_processing/tem"
}

public_s3_base_url <- function(bucket) {
  bucket <- bucket %>% str_remove("/+$")
  if (str_detect(bucket, "^https?://")) {
    return(bucket)
  }

  ## the mc alias (e.g. s3test) is dropped, the rest of the path is served by the public endpoint
  bucket_path <- bucket %>% str_remove("^[^/]+/?")
  if (bucket_path == "") {
    "https://schwab.buckets.embl.de"
  } else {
    file.path("https://schwab.buckets.embl.de", bucket_path)
  }
}

s3_public_base <- public_s3_base_url(s3_bucket)

expected_datasets <- opt$expected_datasets
if (is.null(expected_datasets) || is.na(expected_datasets) || !file.exists(expected_datasets)) {
  stop("--expected_datasets must point to the raw-derived all_datasets.tsv file.")
}

expected_datasets_table <- read_tsv(
  expected_datasets,
  col_types=cols(.default=col_character())
)

expected_omezarr_names <- expected_datasets_table %>%
  select(any_of(c("omezarr_name", "s3_zarr_name"))) %>%
  unlist(use.names=FALSE) %>%
  as.character() %>%
  discard(is.na) %>%
  discard(~.x == "") %>%
  unique()

if ("filename" %in% names(expected_datasets_table)) {
  expected_omezarr_names <- c(
    expected_omezarr_names,
    paste0(expected_datasets_table$filename, "_correctionblend_gradientcorrected.zarr")
  ) %>%
    discard(is.na) %>%
    discard(~.x == "") %>%
    unique()
}

if (length(expected_omezarr_names) == 0) {
  stop("No expected OME-Zarr or Zarr dataset names were found in the raw-derived dataset inventory.")
}

parse_mc_ls_path <- function(line) {
  parsed <- str_match(line, "^\\[.*?\\]\\s+\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2]
  ifelse(
    is.na(parsed),
    str_match(line, "^\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2],
    parsed
  )
}

col_table <- 
  read_lines(opt$all_s3) %>%
  
  
  
  as_tibble() %>%
  filter(!str_detect(value, "coarse_mask")) %>%
  mutate(
    s3_raw=parse_mc_ls_path(value),
    source_name=str_remove(s3_raw, "/$") %>% basename(),
    site=str_extract(source_name, "ATH|BAR|KRI|TAL|NAP|BIL|POR"),
    cell_id=str_extract(source_name, "c0\\d+"),
    size_frac=str_extract(source_name, "\\d+to\\d+"),
    sampling_time=str_extract(source_name, "_(AM|PM|MID|TARA)_") %>% str_remove_all("_"),
    uri=file.path(s3_public_base, s3_raw),
    name=paste0(str_split(source_name, cell_id) %>% map_chr(.,1), cell_id),
    view=site,
    grid=site,
    exclusive=TRUE
    #zarr_root=str_extract(s3_raw, ".*?(_coarse_mask\\.ome\\.zarr|\\.ome\\.zarr|\\.zarr)(?=/|$)")
  ) %>%
  filter(source_name %in% expected_omezarr_names) %>%
  select(
    uri, name, view, grid,
    site, cell_id, size_frac, sampling_time, source_name, exclusive
  )
  #%>%
  # filter(
  #   !is.na(zarr_root),
  #   !str_detect(zarr_root, "^0/")
  # ) %>%
  # distinct(zarr_root) %>%
  # mutate(
  #   s3_raw=zarr_root,
  #   object_name=basename(s3_raw),
  #   source_name=case_when(
  #     str_detect(object_name, "_coarse_mask\\.ome\\.zarr$") ~ str_remove(object_name, "_coarse_mask\\.ome\\.zarr$"),
  #     str_detect(object_name, "\\.ome\\.zarr$") ~ str_remove(object_name, "\\.ome\\.zarr$"),
  #     str_detect(object_name, "_omezarr$") ~ str_remove(object_name, "_omezarr$"),
  #     str_detect(object_name, "\\.zarr$") ~ str_remove(object_name, "\\.zarr$"),
  #     TRUE ~ object_name
  #   ),
  #   source_name=str_remove(source_name, "_correctionblend$"),
  #   source_name=str_remove(source_name, "_gradientcorrected$"),
  #   is_mask=str_detect(object_name, "_coarse_mask\\.ome\\.zarr$"),
  #   uri=file.path("https://s3.embl.de/temscreen", s3_raw),
  #   name=if_else(is_mask, paste0(source_name, " coarse mask"), source_name),
  #   type="intensities",
  #   view=source_name,
  #   display=if_else(is_mask, "coarse cell mask", "TEM image"),
  #   blend=if_else(is_mask, "alpha", "sum"),
  #   color=if_else(is_mask, "magenta", "white"),
  #   contrast_limits=if_else(is_mask, "(0,1)", ""),
  #   mask_kind=if_else(is_mask, "coarse_binary_mask", ""),
  #   format="OmeZarr",
  #   site=str_extract(source_name, "ATH|BAR|KRI|TAL|NAP|BIL|POR"),
  #   cell_id=str_extract(source_name, "c0\\d+"),
  #   size_frac=str_extract(source_name, "\\d+to\\d+"),
  #   sampling_time=str_extract(source_name, "_(AM|PM|MID|TARA)_") %>% str_remove_all("_"),
  #   group=site
  # ) %>%
  # distinct(uri, .keep_all=TRUE) %>%
  # arrange(source_name, is_mask) %>%
  # mutate(
  #   grid_index=match(source_name, unique(source_name)) - 1,
  #   grid_cols=ceiling(sqrt(n_distinct(source_name))),
  #   grid_x=grid_index %% grid_cols,
  #   grid_y=floor(grid_index / grid_cols),
  #   grid="all_images",
  #   grid_position=paste0("(", grid_x, ",", grid_y, ")")
  # ) %>%
  # arrange(source_name, is_mask) %>%


if (nrow(col_table) == 0) {
  stop("No top-level OME-Zarr datasets found in S3 listing.")
}

image_stats_files <- list.files(
  image_stats_dir,
  pattern="_image_stats\\.tsv$",
  full.names=TRUE
)

if (length(image_stats_files) > 0) {
  image_stats <- image_stats_files %>%
    map_dfr(~read_tsv(.x, col_types=cols(.default = col_character()))) %>%
    distinct(name, .keep_all=TRUE)

  col_table <- col_table %>%
    left_join(
      image_stats %>%
        select(name, min_gray, max_gray, contrast_limits),
      by="name"
    )
} else {
  warning("No *_image_stats.tsv files found; collection table was written without contrast limits.")
}

read_image_log <- function() {
  if (sheet_mode == "google") {
    library(googlesheets4)
    json_key <- opt$google_key
    if (is.null(json_key) || is.na(json_key)) {
      stop("--google_key is required when --sheet_mode google")
    }
    gs4_auth(path=json_key)
    read_sheet(image_log_url, sheet=image_log_sheet, col_types="c")
  } else if (file.exists(local_image_log)) {
    read_tsv(local_image_log, col_types=cols(.default = col_character()))
  } else {
    tibble(shortname=character(), site=character())
  }
}

image_log <- read_image_log()

if (nrow(image_log) > 0 && all(c("shortname", "site") %in% names(image_log))) {
  col_table <- col_table %>%
    left_join(
      image_log %>%
        distinct(shortname, site, .keep_all=TRUE),
      by=c("name"="shortname", "site"="site")
    )
} else {
  warning("Image log has no rows or is missing shortname/site columns; collection table was written without annotations.")
}


if (sheet_mode == "google") {
  library(googlesheets4)
  library(googledrive)

  json_key <- opt$google_key
  if (is.null(json_key) || is.na(json_key)) {
    stop("--google_key is required when --sheet_mode google")
  }

  gs4_auth(path=json_key)
  drive_auth(path = json_key)
  write_sheet(col_table, ss = collection_table, sheet=collection_table_sheet)
} else {
  write_tsv(col_table, file=local_collection_table)
}

write_tsv(tibble(done="done"), file="done.tsv")


