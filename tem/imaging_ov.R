library(tidyverse)
library(getopt)
library(cowplot)

spec <- matrix(c(
  # long option                  short  arg  type
  "rawdir",                   "r",   1,   "character",
  "pngdir",        "p",   1,   "character",
  "dryrun", "d",   1,   "character",
  "dryrun_n", "n", 1, "integer",
  "script_dir", "s", 1, "character",
  "sheet_mode", "m", 1, "character",
  "sheet_url", "u", 1, "character",
  "google_key", "k", 1, "character",
  "local_log", "l", 1, "character",
  "existing_s3", "e", 1, "character"
),
ncol = 4,
byrow = TRUE)
opt <- getopt(spec)

arg_file <- commandArgs(FALSE) %>%
  .[str_detect(., "^--file=")] %>%
  str_remove("^--file=") %>%
  .[1]

script_dir <- opt$script_dir
if (is.null(script_dir) || is.na(script_dir)) {
  script_dir <- dirname(normalizePath(arg_file))
}

sheet_mode <- opt$sheet_mode
if (is.null(sheet_mode) || is.na(sheet_mode)) {
  sheet_mode <- "local"
}

trec_tem_googledoc <- opt$sheet_url
if (is.null(trec_tem_googledoc) || is.na(trec_tem_googledoc)) {
  trec_tem_googledoc <- "https://docs.google.com/spreadsheets/d/143uVeeJ72SQE5eK01lzWYCEiT7pJUF3lX7hJl3R9s9I/edit?gid=258669282#gid=258669282"
}

local_log <- opt$local_log
if (is.null(local_log) || is.na(local_log)) {
  local_log <- file.path(dirname(opt$pngdir), "image_log_local.tsv")
}

# opt <- tibble(
#    rawdir="/g/schwab/tem_screen/raw",
#    pngdir="/g/schwab/tem_screen/pngs",
#   #rawdir="/scratch/rheinnec/tem_screen/raw",
#   #pngdir="/scratch/rheinnec/tem_screen/pngs",
#   dryrun="TRUE"
# )

#print(opt$dryrun)

raw_dir <- opt$rawdir
png_dir <- opt$pngdir

existing_s3 <- opt$existing_s3
if (is.null(existing_s3) || is.na(existing_s3)) {
  existing_s3 <- NULL
}

dryrun_n <- opt$dryrun_n
if (is.null(dryrun_n) || is.na(dryrun_n)) {
  dryrun_n <- 10L
}

parse_mc_ls_path <- function(line) {
  parsed <- str_match(line, "^\\[.*?\\]\\s+\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2]
  ifelse(
    is.na(parsed),
    str_match(line, "^\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2],
    parsed
  )
}

existing_s3_paths <- character()
existing_s3_names <- character()
if (!is.null(existing_s3) && file.exists(existing_s3)) {
  existing_s3_paths <- read_lines(existing_s3) %>%
    parse_mc_ls_path() %>%
    str_remove("/$") %>%
    discard(is.na)

  existing_s3_names <- existing_s3_paths %>%
    basename() %>%
    unique()
}

has_s3_omezarr <- function(filename) {
  if (length(existing_s3_names) == 0) {
    return(FALSE)
  }

  expected_names <- c(
    paste0(filename, "_omezarr"),
    paste0(filename, "_correctionblend_gradientcorrected.zarr"),
    paste0(filename, ".ome.zarr"),
    paste0(filename, ".zarr")
  )

  any(expected_names %in% existing_s3_names)
}

all_files_raw <- 
  tibble(file=list.files(raw_dir, pattern="c\\d+.mrc$", recursive = T, full.names=T) %>%
  .[which(!str_detect(.,"canc_"))]) %>%
 # .[63,] %>%
  # filter(
  #   ## duplicated with id c015 from same block
  #   !2474646_76273
  # ) %>%
  #rowwise() %>%
  mutate(
    mdoc_file=str_replace(file, "\\.mrc$", ".mrc.mdoc"),
    site=str_extract(file, "ATH|BAR|KRI|TAL|NAP|BIL|POR"),
    cell_id=str_extract(file, "c0\\d+\\.mrc$") %>% str_remove("\\.mrc"),
    
    shortname=str_extract(basename(dirname(file)), "^.*Cut\\d+") %>% paste(cell_id, sep="_"),
     
    filename=str_split(basename(dirname(file)), "Cut\\d*_") %>% map_chr(.,2) %>% paste(shortname,., sep="_"),
    justblend_file=file.path(png_dir, site, paste0(filename, "_blend.png")),
    correctionblend_file=file.path(png_dir, site, paste0(filename, "_correctionblend.png")),
    png_export_file=file.path(png_dir, site, paste0(filename, "_correctionblend_gradientcorrected.png")),
    omezarr_name=paste0(filename, "_omezarr"),
    s3_zarr_name=paste0(filename, "_correctionblend_gradientcorrected.zarr"),
    
    #grid=
  ) %>%
  rowwise() %>%
  mutate(
    filesize=file.info(file)$size,
    req_mem=min(max(16, round(20*filesize/10^9)), 128),
    s3_omezarr_present=has_s3_omezarr(filename),
    png_export_present=file.exists(png_export_file),
    needs_s3=!s3_omezarr_present,
    needs_png=!png_export_present,
    ## only the final OME-Zarr on S3 decides, the png is kept for information only
    needs_processing=needs_s3
  )

all_files <- all_files_raw %>%
  ungroup() %>%
  select(
    filename, file, mdoc_file, shortname, req_mem,
    justblend_file, correctionblend_file, png_export_file, omezarr_name, s3_zarr_name,
    s3_omezarr_present, png_export_present, needs_s3, needs_png, needs_processing,
    filesize
  ) #%>%

if(as.logical(opt$dryrun)){
  to_run <- all_files %>%
    filter(needs_processing) %>%
    head(dryrun_n)
} else {
  to_run <- all_files %>%
    filter(needs_processing)
}


## make statistics figure

source(file.path(script_dir, "count_stats.R"))

comb_plot <- make_main_statistic_of_sample_number(raw_dir)
outdir <- getwd()
plot_name <- "TEM_screen_image_count"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=4.5, height=4.5)
comb_plot
dev.off()




png(
  filename = file.path(outdir, paste0(plot_name, ".png")),
  width = 4.5,
  height = 4.5,
  units = "in",
  res = 300   # good default for publication-quality
)

comb_plot
dev.off()






# 
# 
# cache_dir <- "/g/schwab/marco/repos/trec_biobank/tem/gargle_cache"
# 
# 
# 
# options(gargle_oauth_cache = cache_dir,
#         gargle_oauth_email = "marco.rheinnecker@embl.de")
# 
# print(gargle::gargle_oauth_cache())
# print(list.files(cache_dir, recursive = TRUE))
# 
# 
# 
# googledrive::drive_auth(
#   email = "marco.rheinnecker@embl.de",
#   cache = cache_dir,
#   scopes = "https://www.googleapis.com/auth/drive"
# )
# 
# folder <- drive_get(as_id("https://drive.google.com/drive/folders/11T1ozEQ66wFDgfWjJnCW3jpfHTHwpJsW"))
# 
# 
# 
# drive_upload(
#   media = plot_name,
#   path  = folder,
#   name  = plot_name,
#   overwrite=TRUE
# )
# 

if(as.logical(opt$dryrun)){
  sheet_name <- "image_log_test"
} else {
  sheet_name <- "image_log"
}



if (sheet_mode == "google") {
  library(googlesheets4)
  library(googledrive)

  json_key <- opt$google_key
  if (is.null(json_key) || is.na(json_key)) {
    json_key <- file.path(script_dir, "trec-tem-screen-e98a2e03f58b.json")
  }

  gs4_auth(path=json_key)
  drive_auth(path = json_key)
  df_trec_tem_current_state <- read_sheet(trec_tem_googledoc, sheet=sheet_name, col_types="c")
} else {
  if (file.exists(local_log)) {
    df_trec_tem_current_state <- read_tsv(local_log, col_types=cols(.default = col_character()))
  } else {
    df_trec_tem_current_state <- tibble(shortname=character(), site=character())
  }
}


if((nrow(df_trec_tem_current_state)>nrow(all_files_raw))&!as.logical(opt$dryrun)){
  
  stop("existing image_log sheet has more entries than present raw files. Proceeding could mean annotations are lost!! - is this suppossed to be a dryrun? Set --dryrun \"TRUE\"")
  
}


timestamp <- Sys.time() %>% as.character() %>% str_replace_all(" |:|\\.", "_")

write_tsv(df_trec_tem_current_state, file=paste0("manually_filled_log_",timestamp,".tsv"))


new <- all_files_raw %>%
  select(shortname, site) %>%
  left_join(df_trec_tem_current_state, by=c("shortname", "site"))

#, "site"

if (sheet_mode == "google") {
  write_sheet(new, ss = trec_tem_googledoc, sheet=sheet_name)
} else {
  write_tsv(new, file=local_log)
}

write_csv(to_run, file="images_to_process.csv")

write_tsv(all_files, file="all_datasets.tsv")






