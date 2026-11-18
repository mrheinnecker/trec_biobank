library(tidyverse)
library(getopt)

spec <- matrix(c(
  "input_table", "i", 1, "character",
  "sheet_mode", "m", 1, "character",
  "sheet_url", "u", 1, "character",
  "sheet_name", "s", 1, "character",
  "google_key", "k", 1, "character",
  "outdir", "o", 1, "character",
  "dryrun", "d", 1, "character",
  "dryrun_n", "n", 1, "integer",
  "dryrun_images", NA, 1, "character",
  "existing_s3", "e", 1, "character",
  "default_x_scale", NA, 1, "character",
  "default_y_scale", NA, 1, "character",
  "default_z_scale", NA, 1, "character",
  "scale_unit", NA, 1, "character"
), ncol=4, byrow=TRUE)
opt <- getopt(spec)

# Local debug inputs:
# Uncomment this block when running the script manually, and keep it commented
# when Nextflow is passing the real command-line options.
#
# opt <- list(
#   input_table = "C:/projects/exm_screen/exm_images.tsv",
#   sheet_mode = "google",
#   sheet_url = "https://docs.google.com/spreadsheets/d/REPLACE_WITH_EXM_INPUT_SHEET_ID/edit?gid=0#gid=0",
#   sheet_name = "exm_input_table",
#   google_key = "/g/schwab/marco/repos/trec_biobank/trec-tem-screen-e98a2e03f58b.json",
#   outdir = "/scratch/rheinnec/central_data_processing/exm_screen",
#   dryrun = "FALSE",
#   dryrun_n = 2L,
#   existing_s3 = "C:/projects/exm_screen/existing_s3_entries.txt",
#   default_x_scale = "",
#   default_y_scale = "",
#   default_z_scale = "",
#   scale_unit = "nm"
# )

value_or_default <- function(value, default) {
  if (is.null(value) || is.na(value) || value == "") default else value
}

sheet_mode <- value_or_default(opt$sheet_mode, "local")
sheet_name <- value_or_default(opt$sheet_name, "")
dryrun <- value_or_default(opt$dryrun, "FALSE")
dryrun_n <- opt$dryrun_n
if (is.null(dryrun_n) || is.na(dryrun_n)) {
  dryrun_n <- 2L
}
outdir <- value_or_default(opt$outdir, "exm_processed")
scale_unit <- value_or_default(opt$scale_unit, "nm")

parse_mc_ls_path <- function(line) {
  parsed <- str_match(line, "^\\[.*?\\]\\s+\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2]
  ifelse(
    is.na(parsed),
    ## plain names (listing of the mounted S3 share) are used as they are
    coalesce(str_match(line, "^\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2], line),
    parsed
  )
}

existing_s3_names <- character()
if (!is.null(opt$existing_s3) && !is.na(opt$existing_s3) && file.exists(opt$existing_s3)) {
  existing_s3_paths <- read_lines(opt$existing_s3) %>%
    parse_mc_ls_path() %>%
    str_remove("/$") %>%
    discard(is.na)

  existing_s3_names <- existing_s3_paths %>%
    str_match("^([^/]+?)(?:\\.ome\\.zarr|\\.zarr)$") %>%
    .[, 2] %>%
    discard(is.na) %>%
    unique()
}

sanitize_name <- function(x) {
  tools::file_path_sans_ext(basename(x)) %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

read_local_table <- function(path) {
  if (is.null(path) || is.na(path)) {
    stop("--input_table is required when --sheet_mode local")
  }
  if (str_detect(path, "\\.[Cc][Ss][Vv]$")) {
    read_csv(path, col_types=cols(.default=col_character()))
  } else {
    read_tsv(path, col_types=cols(.default=col_character()))
  }
}

read_google_table <- function(url, sheet) {
  if (is.null(url) || is.na(url) || url == "") {
    stop("--sheet_url is required when --sheet_mode google")
  }
  json_key <- opt$google_key
  if (is.null(json_key) || is.na(json_key)) {
    stop("--google_key is required when --sheet_mode google")
  }

  library(googlesheets4)
  gs4_auth(path=json_key)
  if (is.null(sheet) || is.na(sheet) || sheet == "") {
    read_sheet(url, col_types="c")
  } else {
    read_sheet(url, sheet=sheet, col_types="c")
  }
}

first_existing_column <- function(data, candidates) {
  matches <- intersect(candidates, names(data))
  if (length(matches) == 0) NULL else matches[[1]]
}

coalesce_column_values <- function(data, candidates, default="") {
  column <- first_existing_column(data, candidates)
  if (is.null(column)) {
    rep(default, nrow(data))
  } else {
    coalesce(data[[column]], default)
  }
}

images_raw <- if (sheet_mode == "google") {
  read_google_table(opt$sheet_url, sheet_name)
} else {
  read_local_table(opt$input_table)
}

images <- images_raw %>%
  mutate(across(everything(), ~na_if(.x, "")))

path_column <- first_existing_column(
  images,
  c("FilePath", "raw_path", "file_path", "filepath", "file", "source_path", "path")
)
if (is.null(path_column)) {
  stop("Input table must contain a FilePath, raw_path, file_path, filepath, file, source_path, or path column")
}

if (!"convert" %in% names(images)) {
  images[["convert"]] <- "1"
}

images[["exm_x_scale"]] <- coalesce_column_values(
  images,
  c("x_scale", "x_pixel_size_nm", "pixel_size_x_nm", "physical_size_x_nm"),
  value_or_default(opt$default_x_scale, "")
)
images[["exm_y_scale"]] <- coalesce_column_values(
  images,
  c("y_scale", "y_pixel_size_nm", "pixel_size_y_nm", "physical_size_y_nm"),
  value_or_default(opt$default_y_scale, "")
)
images[["exm_z_scale"]] <- coalesce_column_values(
  images,
  c("z_scale", "z_distance_nm", "z_spacing_nm", "pixel_size_z_nm", "physical_size_z_nm"),
  value_or_default(opt$default_z_scale, "")
)
images[["exm_scale_unit"]] <- coalesce_column_values(
  images,
  c("scale_unit", "pixel_size_unit", "physical_size_unit"),
  scale_unit
)

## scale overrides must be numbers, anything else (e.g. TRUE from another sheet column) is ignored
## and the pixel size from the .nd2 file is used
numeric_or_empty <- function(value, column, image) {
  if (is.na(value) || value == "") {
    return("")
  }
  if (is.na(suppressWarnings(as.numeric(value)))) {
    warning("Ignoring non-numeric ", column, " '", value, "' for ", image)
    return("")
  }
  value
}

## a row can point to a single .nd2 file or to a folder with .nd2 files
expand_nd2_paths <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern="\\.nd2$", ignore.case=TRUE, full.names=TRUE)
    if (length(files) == 0) {
      warning("No .nd2 files found in ", path)
    }
    files
  } else {
    path
  }
}

## the Nikon software writes one file per channel:
##   Img0098_ChannelCy5 (Single EM)_Seq0068.nd2
## -> all files with the same folder + image prefix (Img0098) form one dataset with one channel each
## tolerated variants: "Img0126 _Channel...", "Img0028__Channel...", missing "_Seq0000", copies ending in " (1)"
channel_file_pattern <- "^(.*?)[\\s_]*_Channel(.+?)(?:_Seq\\d+)?(?:\\s*\\(\\d+\\))?\\s*$"

files <- images %>%
  mutate(
    raw_file=str_remove(.data[[path_column]], "/+$")
  ) %>%
  filter(!is.na(raw_file), raw_file != "") %>%
  mutate(raw_file=map(raw_file, expand_nd2_paths)) %>%
  unnest_longer(raw_file) %>%
  distinct(raw_file, .keep_all=TRUE) %>%
  mutate(
    file_stem=tools::file_path_sans_ext(basename(raw_file)),
    image_prefix=str_trim(coalesce(str_match(file_stem, channel_file_pattern)[, 2], file_stem)),
    channel_name=str_trim(str_match(file_stem, channel_file_pattern)[, 3]),
    ## "Cy5 (Single EM)" -> "Cy5", to find the same channel saved twice
    channel_key=str_trim(str_remove_all(channel_name, "\\(.*?\\)")),
    image_dir=dirname(raw_file),
    ## image numbers (Img0098) restart per session -> the folder name is part of the dataset name
    default_filename=paste(basename(image_dir), image_prefix, sep="_")
  )

all_images <- files %>%
  group_by(image_dir, image_prefix) %>%
  summarise(
    raw_path=paste(sort(raw_file), collapse=";"),
    n_files=n(),
    channel_names=paste(sort(na.omit(channel_name)), collapse=";"),
    ## dataset name: dataset_name column, else shortname (e.g. BAR_Img0093), else <folder>_<image prefix>
    ## (a filename column holds the raw file names and is not used)
    filename=first(na.omit(c(
      if ("dataset_name" %in% names(files)) dataset_name,
      if ("shortname" %in% names(files)) shortname,
      default_filename
    ))),
    shortname=if ("shortname" %in% names(files)) first(na.omit(c(shortname, image_prefix))) else first(image_prefix),
    req_mem=if ("req_mem" %in% names(files)) first(na.omit(c(req_mem, "32"))) else "32",
    x_scale=numeric_or_empty(first(exm_x_scale), "x_scale", first(image_prefix)),
    y_scale=numeric_or_empty(first(exm_y_scale), "y_scale", first(image_prefix)),
    z_scale=numeric_or_empty(first(exm_z_scale), "z_scale", first(image_prefix)),
    scale_unit=first(exm_scale_unit),
    convert=first(convert),
    .groups="drop"
  ) %>%
  mutate(
    filename=sanitize_name(paste0(filename, ".nd2")),
    output_path=file.path(outdir, filename, paste0(filename, ".ome.zarr")),
    s3_omezarr_present=filename %in% existing_s3_names,
    convert_selected=!tolower(coalesce(convert, "1")) %in% c("0", "false", "no", "n"),
    needs_processing=convert_selected & !s3_omezarr_present
  ) %>%
  select(filename, shortname, raw_path, output_path, req_mem, x_scale, y_scale, z_scale, scale_unit, everything())

## the same channel in two files (e.g. a copy "..._Seq0016 (1).nd2") would end up as an extra channel
double_channels <- files %>%
  filter(!is.na(channel_key)) %>%
  group_by(image_dir, image_prefix, channel_key) %>%
  filter(n() > 1) %>%
  ungroup()
if (nrow(double_channels) > 0) {
  stop(
    "The same channel is listed in more than one file of an image; keep one of them in the input table:\n  ",
    paste(double_channels$raw_file, collapse="\n  ")
  )
}

## all channel rows of one image must carry the same name, otherwise the first file's value would be used silently
name_columns <- intersect(c("dataset_name", "shortname"), names(files))
inconsistent <- files %>%
  group_by(image_dir, image_prefix) %>%
  summarise(across(all_of(name_columns), ~n_distinct(na.omit(.x))), .groups="drop") %>%
  filter(if_any(all_of(name_columns), ~.x > 1))
if (nrow(inconsistent) > 0) {
  stop(
    "The channel rows of an image have different ", paste(name_columns, collapse="/"), " values: ",
    paste(file.path(inconsistent$image_dir, inconsistent$image_prefix), collapse=", "),
    ". Use one name per image (all 4 channel rows the same), e.g. shortname BAR_Img0093; leave dataset_name out."
  )
}

if (any(duplicated(all_images$filename))) {
  duplicates <- all_images %>%
    filter(filename %in% filename[duplicated(filename)]) %>%
    mutate(entry=paste0(filename, " <- ", image_dir, "/", image_prefix))
  stop(
    "The same dataset name is used for different images (Nikon image numbers restart per session):\n  ",
    paste(duplicates$entry, collapse="\n  "),
    "\nMake the shortname unique, e.g. add the date or sample: BAR_20240320_Img0093."
  )
}

to_run <- all_images %>%
  filter(needs_processing)
if (as.logical(dryrun)) {
  ## dryrun_images (comma-separated dataset names, e.g. BAR_Img0098,BAR_Img0099) are taken first,
  ## the rest of dryrun_n is filled in table order
  dryrun_images <- value_or_default(opt$dryrun_images, "")
  dryrun_images <- if (tolower(dryrun_images) == "true") character() else str_trim(str_split_1(dryrun_images, ","))
  dryrun_images <- dryrun_images[dryrun_images != ""]
  not_selectable <- setdiff(dryrun_images, to_run$filename)
  if (length(not_selectable) > 0) {
    warning(
      "dryrun_images not in the images to process (unknown name or already published): ",
      paste(not_selectable, collapse=", ")
    )
  }
  to_run <- to_run %>%
    arrange(!filename %in% dryrun_images) %>%
    head(dryrun_n)
}

write_csv(to_run, file="images_to_process.csv")
write_tsv(all_images, file="all_datasets.tsv")
