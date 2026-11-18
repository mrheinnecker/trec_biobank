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
#   input_table = "C:/projects/cryo_screen/cryo_images.tsv",
#   sheet_mode = "google",
#   sheet_url = "https://docs.google.com/spreadsheets/d/1ePRpa56mmMvCeRTLXmwOywOLy5_I3AFrxJepSUYGR1s/edit?gid=1442254503#gid=1442254503",
#   sheet_name = "cryo_lm",
#   google_key = "/g/schwab/marco/repos/trec_biobank/trec-tem-screen-e98a2e03f58b.json",
#   outdir = "/scratch/rheinnec/central_data_processing/cryo_screen",
#   dryrun = "FALSE",
#   dryrun_n = 2L,
#   existing_s3 = "C:/projects/cryo_screen/existing_s3_entries.txt",
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
outdir <- value_or_default(opt$outdir, "cryo_processed")
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

images <- images_raw #%>%
  # mutate(across(everything(), ~na_if(.x, ""))) %>%
  # filter(str_detect(`File Name`, "_st_3D")) %>% 
  # filter(str_detect(FilePath, "CryoLM")) 

path_column <- first_existing_column(
  images,
  c("FilePath","raw_path", "file_path", "filepath", "file", "source_path", "path")
)
if (is.null(path_column)) {
  stop("Input table must contain a raw_path, file_path, filepath, file, source_path, or path column")
}

if (!"convert" %in% names(images)) {
  images[["convert"]] <- "1"
}

images[["cryo_x_scale"]] <- coalesce_column_values(
  images,
  c("x_scale", "x_pixel_size_nm", "pixel_size_x_nm", "physical_size_x_nm"),
  value_or_default(opt$default_x_scale, "")
)
images[["cryo_y_scale"]] <- coalesce_column_values(
  images,
  c("y_scale", "y_pixel_size_nm", "pixel_size_y_nm", "physical_size_y_nm"),
  value_or_default(opt$default_y_scale, "")
)
images[["cryo_z_scale"]] <- coalesce_column_values(
  images,
  c("z_scale", "z_distance_nm", "z_spacing_nm", "pixel_size_z_nm", "physical_size_z_nm"),
  value_or_default(opt$default_z_scale, "")
)
images[["cryo_scale_unit"]] <- coalesce_column_values(
  images,
  c("scale_unit", "pixel_size_unit", "physical_size_unit"),
  scale_unit
)

all_images <- images %>%
  mutate(
    raw_path=.data[[path_column]],
    raw_path=str_remove(raw_path, "/+$")
  ) %>%
  filter(!is.na(raw_path), raw_path != "") %>%
  mutate(
    filename=if ("filename" %in% names(images)) coalesce(filename, sanitize_name(raw_path)) else sanitize_name(raw_path),
    filename=sanitize_name(filename),
    shortname=if ("shortname" %in% names(images)) coalesce(shortname, filename) else filename,
    output_path=file.path(outdir, filename, paste0(filename, ".ome.zarr")),
    req_mem=if ("req_mem" %in% names(images)) coalesce(req_mem, "32") else "32",
    x_scale=cryo_x_scale,
    y_scale=cryo_y_scale,
    z_scale=cryo_z_scale,
    scale_unit=cryo_scale_unit,
    s3_omezarr_present=filename %in% existing_s3_names,
    convert_selected=!tolower(coalesce(convert, "1")) %in% c("0", "false", "no", "n"),
    needs_processing=convert_selected & !s3_omezarr_present
  ) %>%
  distinct(raw_path, .keep_all=TRUE) %>%
  select(filename, shortname, raw_path, output_path, req_mem, x_scale, y_scale, z_scale, scale_unit, everything())

to_run <- all_images %>%
  filter(needs_processing)
if (as.logical(dryrun)) {
  to_run <- head(to_run, dryrun_n)
}

write_csv(to_run, file="images_to_process.csv")
write_tsv(all_images, file="all_datasets.tsv")
