

library(tidyverse)

ov_dir <- "C:/projects/tem_classification/ov"

## set TRUE to export the sheets again (needs Google auth), otherwise the exported TSVs are used
refresh_from_sheets <- TRUE

if (refresh_from_sheets) {
  sampling_events <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1jcpyMkSR4npSxST3D5cFzkAIi9UPmwbPvzAdr2ws55U/edit?gid=193823038#gid=193823038", sheet="Event", col_types = "c", skip=1)

  
  cryo_sheet <- "https://docs.google.com/spreadsheets/d/1Ln2YVfxMtQNc-eSh69RA0tpc5UhsMNFXizfCVAR4iZI/edit?gid=0#gid=0"
  cryo_pf <- googlesheets4::read_sheet(cryo_sheet, sheet="pf", col_types="c")
  cryo_wm <- googlesheets4::read_sheet(cryo_sheet, sheet="wm", col_types="c")
  hpf <- googlesheets4::read_sheet(cryo_sheet, sheet="hpf", col_types="c")
  
  write_tsv(sampling_events, file=file.path(ov_dir, "sampling_events.tsv"))
  write_tsv(hpf, file=file.path(ov_dir, "hpf.tsv"))
  write_tsv(cryo_pf, file=file.path(ov_dir, "cryo_pf.tsv"))
  write_tsv(cryo_wm, file=file.path(ov_dir, "cryo_wm.tsv"))
  
  
 # hpf <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1jcpyMkSR4npSxST3D5cFzkAIi9UPmwbPvzAdr2ws55U/edit?gid=193823038#gid=193823038", sheet="HPF", col_types = "c")

  write_tsv(sampling_events, file=file.path(ov_dir, "sampling_events.tsv"))
  write_tsv(hpf, file=file.path(ov_dir, "hpf.tsv"))


  ## metabarcoding samples (same sheet as Event / HPF)
  metab <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1jcpyMkSR4npSxST3D5cFzkAIi9UPmwbPvzAdr2ws55U/edit?gid=804034428#gid=804034428", sheet="MetaB", col_types="c")
  write_tsv(metab, file=file.path(ov_dir, "metab.tsv"))

  ## chemically fixed samples (same sheet as Event / HPF)
  chemfix <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1jcpyMkSR4npSxST3D5cFzkAIi9UPmwbPvzAdr2ws55U/edit?gid=193823038#gid=193823038", sheet="Chem fixation", col_types="c")
  write_tsv(chemfix, file=file.path(ov_dir, "chemfix.tsv"))
}


## ---- imaging: MoBIE collection tables written by the workflows ---------------------------------
## one tab per workflow in the collection table sheet; saved as <ov_dir>/collection_tables/<tab>.tsv
refresh_collection_tables <- TRUE

collection_table_sheet <- "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit"
collection_table_tabs <- c(
  "tem_collection_table",
  "sem_collection_table",
  "hitt_collection_table",
  "cryo_collection_table",
  "cryo_hpf_collection_table",
  "plastic_collection_table",
  "exm_collection_table"
)
collection_table_dir <- file.path(ov_dir, "collection_tables")

if (refresh_collection_tables) {
  dir.create(collection_table_dir, showWarnings=FALSE, recursive=TRUE)
  ## tabs in the sheet, to report missing ones instead of failing
  available_tabs <- googlesheets4::sheet_names(collection_table_sheet)
  message("Tabs in the collection table sheet: ", paste(available_tabs, collapse=", "))

  walk(collection_table_tabs, function(tab) {
    if (!tab %in% available_tabs) {
      warning("Collection table tab not found: ", tab)
      return(invisible(NULL))
    }
    tab_data <- googlesheets4::read_sheet(collection_table_sheet, sheet=tab, col_types="c")
    write_tsv(tab_data, file=file.path(collection_table_dir, paste0(tab, ".tsv")))
    message(tab, ": ", nrow(tab_data), " rows")
  })
}

collection_tables <- list.files(collection_table_dir, pattern="\\.tsv$", full.names=TRUE) %>%
  set_names(str_remove(basename(.), "\\.tsv$")) %>%
  map(~read_tsv(.x, col_types=cols(.default=col_character())))


## ---- EXM annotation sheet: all tabs, saved as <ov_dir>/exm_annotations/<tab>.tsv ------------------
refresh_exm_annotations <- TRUE

exm_annotation_sheet <- "https://docs.google.com/spreadsheets/d/1PZzLraZx53w_zKwpISj_r7riUibX6iGg9Z2Li_Gc3rw/edit"
exm_annotation_dir <- file.path(ov_dir, "exm_annotations")

if (refresh_exm_annotations) {
  dir.create(exm_annotation_dir, showWarnings=FALSE, recursive=TRUE)
  exm_tabs <- googlesheets4::sheet_names(exm_annotation_sheet)
  message("Tabs in the EXM annotation sheet: ", paste(exm_tabs, collapse=", "))
  walk(exm_tabs, function(tab) {
    tab_data <- googlesheets4::read_sheet(exm_annotation_sheet, sheet=tab, col_types="c")
    file_name <- paste0(str_replace_all(tab, "[^A-Za-z0-9_-]+", "_"), ".tsv")
    write_tsv(tab_data, file=file.path(exm_annotation_dir, file_name))
    message(tab, ": ", nrow(tab_data), " rows")
  })
}

sampling_events <- read_tsv(file.path(ov_dir, "sampling_events.tsv"), col_types=cols(.default=col_character()))
hpf <- read_tsv(file.path(ov_dir, "hpf.tsv"), col_types=cols(.default=col_character()))
cryo_pf <- read_tsv(file.path(ov_dir, "cryo_pf.tsv"), col_types=cols(.default=col_character()))
cryo_wm <- read_tsv(file.path(ov_dir, "cryo_wm.tsv"), col_types=cols(.default=col_character()))
metab <- read_tsv(file.path(ov_dir, "metab.tsv"), col_types=cols(.default=col_character()))
chemfix <- read_tsv(file.path(ov_dir, "chemfix.tsv"), col_types=cols(.default=col_character()))


## ---- samples -> sampling events -----------------------------------------------------------------
## no shared ID column: matched by site (sample "Site" = event "Country") and date,
## on days with several events additionally by the sample time of day (AM / MID / PM).
##
## match_quality:
##   date          only one event at the site on that day
##   date_ampm     several events, exactly one in the same half day (MID: 10-14h)
##   untimed_event several events, none in the same half day, exactly one without collection time
##   nearest_time  still ambiguous: the event closest to AM 09:00 / MID 12:00 / PM 15:00 (check these)
##   no_sample_time several events but the sample has no time of day: not assigned, see candidate_events
##   no_event      no event at the site on that day (or no valid date)

events <- sampling_events %>%
  transmute(
    event_id=`Event ID`,
    event_site_id=`Site ID`,
    site=Country,
    date=mdy(`Collection date`),
    event_time=na_if(`Collection time`, "NA"),
    event_hour=as.numeric(str_sub(event_time, 1, 2)) + as.numeric(str_sub(event_time, 4, 5)) / 60,
    latitude=as.numeric(Latitude),
    longitude=as.numeric(Longitude)
  )

target_hour <- c(AM=9, MID=12, PM=15)

in_half_day <- function(hour, time_of_day) {
  case_when(
    is.na(hour) | is.na(time_of_day) ~ NA,
    time_of_day == "AM" ~ hour < 12,
    time_of_day == "PM" ~ hour >= 12,
    time_of_day == "MID" ~ hour >= 10 & hour <= 14,
    TRUE ~ NA
  )
}

## samples: a table with site, date (Date), time_of_day; returns it with the matched event columns
map_to_events <- function(samples, events) {
  ## one row per time slot (site, date, time of day) and candidate event
  slot_event <- samples %>%
    distinct(site, date, time_of_day) %>%
    filter(!is.na(date)) %>%
    inner_join(events, by=c("site", "date"), relationship="many-to-many") %>%
    group_by(site, date, time_of_day) %>%
    mutate(
      n_events=n(),
      same_half=in_half_day(event_hour, time_of_day),
      n_same_half=sum(same_half %in% TRUE),
      n_untimed=sum(is.na(event_hour)),
      time_distance=abs(event_hour - target_hour[time_of_day]),
      candidate_events=paste0(event_site_id, " (", coalesce(str_sub(event_time, 1, 5), "no time"), ")", collapse="; "),
      match_quality=case_when(
        n_events == 1 ~ "date",
        is.na(time_of_day) ~ "no_sample_time",
        n_same_half == 1 ~ "date_ampm",
        n_same_half == 0 & n_untimed == 1 ~ "untimed_event",
        TRUE ~ "nearest_time"
      ),
      keep=case_when(
        match_quality == "date" ~ TRUE,
        match_quality == "date_ampm" ~ same_half %in% TRUE,
        match_quality == "untimed_event" ~ is.na(event_hour),
        match_quality == "no_sample_time" ~ row_number() == 1,
        ## several in the same half day: the closest of those, else the closest overall (untimed last)
        n_same_half > 1 ~ same_half %in% TRUE & time_distance == suppressWarnings(min(time_distance[same_half %in% TRUE])),
        TRUE ~ !is.na(time_distance) & time_distance == suppressWarnings(min(time_distance, na.rm=TRUE))
      )
    ) %>%
    filter(keep) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(across(c(event_id, event_site_id, event_time), ~if_else(match_quality == "no_sample_time", NA, .x))) %>%
    mutate(across(c(latitude, longitude), ~if_else(match_quality == "no_sample_time", NA_real_, .x))) %>%
    select(site, date, time_of_day, event_id, event_site_id, event_time, latitude, longitude, n_events, candidate_events, match_quality)

  samples %>%
    left_join(slot_event, by=c("site", "date", "time_of_day")) %>%
    mutate(match_quality=coalesce(match_quality, "no_event"))
}

## overview of a mapping; the flagged slots should be checked by hand
report_mapping <- function(mapped, label) {
  message("---- ", label)
  print(count(mapped, match_quality))
  mapped %>%
    filter(match_quality %in% c("untimed_event", "nearest_time", "no_sample_time", "no_event")) %>%
    count(site, date_raw, time_of_day, event_time, event_site_id, candidate_events, match_quality) %>%
    print(n=Inf)
  invisible(mapped)
}

sample_base <- function(tab) {
  tab %>%
    mutate(
      barcode=Barcode,
      site=Site,
      date=suppressWarnings(ymd(Date)),
      date_raw=Date,
      time_of_day=na_if(Time, "NA"),
      tara=TARA
    )
}


## ---- HPF ----
hpf_samples <- hpf %>%
  ## the column was renamed in the sheet
  rename(any_of(c(block_name_if_produced="Block name (if produced)"))) %>%
  sample_base() %>%
  transmute(
    barcode, site, date, date_raw, time_of_day, tara,
    parent_barcode=na_if(`Parent HPF barcode (if at all)`, "NA"),
    ## typos in the sheet: "NAP _10to40_..." and a " CHOPPED" suffix
    block_name=na_if(block_name_if_produced, "NA") %>%
      str_remove(" CHOPPED$") %>%
      str_remove_all(" "),
    cryo_name=na_if(`Cryo name (automatically generated FP)`, "NA"),
    used=na_if(Used, "NA"),
    imaged=Imaged,
    post_processed=`Post Processed`,
    freezing_quality=na_if(eye_test_freezing_quality, "NA"),
    amount_of_signal=na_if(amount_of_signal, "NA")
  )

hpf_events <- map_to_events(hpf_samples, events) %>%
  report_mapping("HPF")

write_tsv(hpf_events, file=file.path(ov_dir, "hpf_event_map.tsv"))


## ---- cryo plunge-frozen (pf) and wm ----
cryo_pf_samples <- cryo_pf %>%
  sample_base() %>%
  transmute(
    barcode, site, date, date_raw, time_of_day, tara,
    fraction=Fraction,
    name=`New name`,
    imaged=`Imaged with Linkam Stage`,
    post_processed=`Post Processed`,
    freezing_quality=na_if(Freezing_quality, "NA"),
    sample_quality=na_if(Sample_quality, "NA")
  )

cryo_wm_samples <- cryo_wm %>%
  sample_base() %>%
  transmute(
    barcode, site, date, date_raw, time_of_day, tara,
    fraction=Fraction,
    name=`New name`,
    imaged=`Imaged with Linkam Stage Cryo-WF`,
    comments=na_if(Comments, "NA"),
    post_processed=`Post Processed`,
    freezing_quality=na_if(Freezing_quality, "NA"),
    sample_quality=na_if(Sample_quality, "NA")
  )

cryo_pf_events <- map_to_events(cryo_pf_samples, events) %>%
  report_mapping("cryo PF")
cryo_wm_events <- map_to_events(cryo_wm_samples, events) %>%
  report_mapping("cryo WM")

write_tsv(cryo_pf_events, file=file.path(ov_dir, "cryo_pf_event_map.tsv"))
write_tsv(cryo_wm_events, file=file.path(ov_dir, "cryo_wm_event_map.tsv"))


## ---- metabarcoding ----
metab_samples <- metab %>%
  ## empty rows and note rows at the end of the sheet
  filter(!is.na(Site), Site != "NA") %>%
  rename(Barcode=`Barcode ID`) %>%
  sample_base() %>%
  transmute(
    barcode, site, date, date_raw, time_of_day, tara,
    lab_id=`LabID name`,
    ## "10-40" and "10to40" both occur
    fraction=str_replace(Fraction, "-", "to"),
    replicate=Replicate,
    depth_m=`Depth (m)`,
    sent_to_genoscope=`Sent to Genoscope`,
    notes=na_if(Notes, "NA")
  )

metab_events <- map_to_events(metab_samples, events) %>%
  report_mapping("MetaB")

write_tsv(metab_events, file=file.path(ov_dir, "metab_event_map.tsv"))


## ---- chemical fixation ----
## site as code in this sheet
site_codes <- c(
  ROS="Roscoff", TAL="Tallinn", KRI="Kristineberg", BIL="Bilbao",
  POR="Porto", BAR="Barcelona", NAP="Naples", ATH="Athens"
)

chemfix_samples <- chemfix %>%
  mutate(
    Site=unname(coalesce(site_codes[Site], Site)),
    Barcode=`TREC ID`,
    TARA=TARA_AML,
    ## year typo in the sheet (also in the sample names): Barcelona was sampled in March 2024 only
    ## (all events, HPF and MetaB), the rows dated March 2023 are the same days in 2024
    date_corrected=Site == "Barcelona" & str_detect(Date, "^202303"),
    Date=if_else(date_corrected, str_replace(Date, "^2023", "2024"), Date)
  ) %>%
  sample_base() %>%
  transmute(
    barcode, site, date, date_raw, date_corrected, time_of_day, tara,
    name=`Sample name`,
    ## "20-200" and "20to200" both occur
    fraction=str_replace(`Size fraction`, "-", "to"),
    tube=`Tube number`,
    used=Used
  )

chemfix_events <- map_to_events(chemfix_samples, events) %>%
  report_mapping("Chem fixation")

write_tsv(chemfix_events, file=file.path(ov_dir, "chemfix_event_map.tsv"))


## ---- imaging -> samples / events ------------------------------------------------------------
## images from the collection tables (EXM skipped for now), linked in this order:
##   barcode       9-digit barcode in the image name -> sample of any layer
##   block_name    resin block name in the image name -> HPF sample with that block
##   name_event    site, date and time of day parsed from the name (or columns) -> sampling event only
##   none          nothing found

## all physical samples, one row per barcode
sample_index <- bind_rows(
  hpf_events %>% transmute(barcode, block_name, sample_layer="HPF", event_id),
  cryo_pf_events %>% transmute(barcode, sample_layer="cryo PF", event_id),
  cryo_wm_events %>% transmute(barcode, sample_layer="cryo WM", event_id),
  metab_events %>% transmute(barcode, sample_layer="MetaB", event_id),
  chemfix_events %>% transmute(barcode, sample_layer="Chem fixation", event_id)
) %>%
  ## "-" is used as placeholder in the chem fixation sheet
  filter(!is.na(barcode), !barcode %in% c("NA", "-"))

duplicated_barcodes <- sample_index %>% count(barcode) %>% filter(n > 1)
if (nrow(duplicated_barcodes) > 0) {
  warning(nrow(duplicated_barcodes), " barcodes occur in more than one sample row, the first is used")
}
sample_by_barcode <- sample_index %>%
  distinct(barcode, .keep_all=TRUE) %>%
  select(barcode, sample_layer, sample_event_id=event_id)

sample_by_block <- sample_index %>%
  filter(!is.na(block_name)) %>%
  distinct(block_name, .keep_all=TRUE) %>%
  select(block_name, block_barcode=barcode, block_event_id=event_id)

## SITE_[TARA_]fraction_date_time...: site name, date and time of day for the event fallback
parse_sample_name <- function(x) {
  m <- str_match(x, "(?:^|_)([A-Z]{3})_(?:TARA_)?(?:(\\d+to\\d+|NA)_)?(20\\d{6})_(AM|PM|MID)?")
  tibble(
    name_site=unname(site_codes[m[, 2]]),
    name_date=ymd(m[, 4], quiet=TRUE),
    name_time=m[, 5]
  )
}

link_images <- function(images) {
  images %>%
    bind_cols(parse_sample_name(images$parse_name)) %>%
    left_join(sample_by_barcode, by="barcode") %>%
    left_join(sample_by_block, by="block_name") %>%
    left_join(
      events %>% select(site, date, event_id),
      by=c("name_site"="site", "name_date"="date"),
      relationship="many-to-many"
    ) %>%
    ## several events on that day: decided by the time of day as for the samples
    group_by(image_id) %>%
    mutate(n_name_events=n_distinct(event_id, na.rm=TRUE)) %>%
    ungroup() %>%
    select(-event_id) %>%
    distinct() %>%
    left_join(
      map_to_events(
        distinct(., site=name_site, date=name_date, time_of_day=name_time) %>% filter(!is.na(site)),
        events
      ) %>% select(site, date, time_of_day, name_event_id=event_id),
      by=c("name_site"="site", "name_date"="date", "name_time"="time_of_day")
    ) %>%
    mutate(
      link=case_when(
        !is.na(sample_layer) ~ "barcode",
        !is.na(block_barcode) ~ "block_name",
        !is.na(name_event_id) ~ "name_event",
        TRUE ~ "none"
      ),
      sample_barcode=case_when(link == "barcode" ~ barcode, link == "block_name" ~ block_barcode),
      sample_layer=case_when(link == "barcode" ~ sample_layer, link == "block_name" ~ "HPF"),
      event_id=case_when(
        link == "barcode" ~ sample_event_id,
        link == "block_name" ~ block_event_id,
        link == "name_event" ~ name_event_id
      )
    ) %>%
    select(imaging, image_id, uri, sample_barcode, sample_layer, event_id, link, block_name, fraction, any_of("exm_barcode"), name_site, name_date, name_time)
}

barcode_in <- function(x) str_match(x, "_(\\d{9})_")[, 2]
## one image per uri (the cryo / plastic tables have one row per channel)
one_per_uri <- function(tab) distinct(tab, uri, .keep_all=TRUE)

image_tables <- list()

if ("cryo_collection_table" %in% names(collection_tables)) {
  image_tables$cryo <- collection_tables$cryo_collection_table %>%
    one_per_uri() %>%
    transmute(
      imaging=if_else(str_detect(name, "_PF_"), "cryo LM PF", "cryo LM WM"),
      image_id=str_remove(name, "_c\\d+_[^_]+$"),
      uri, barcode=barcode_in(name), block_name=NA_character_, parse_name=name
    )
}

## cryo confocal screens of HPF carriers (own collection table), linked to the HPF sample by barcode
if ("cryo_hpf_collection_table" %in% names(collection_tables)) {
  image_tables$cryo_hpf <- collection_tables$cryo_hpf_collection_table %>%
    one_per_uri() %>%
    transmute(
      imaging="cryo LM HPF",
      image_id=str_remove(name, "_c\\d+_[^_]+$"),
      uri, barcode=barcode_in(name), block_name=NA_character_, parse_name=name
    )
}

if ("plastic_collection_table" %in% names(collection_tables)) {
  image_tables$plastic <- collection_tables$plastic_collection_table %>%
    one_per_uri() %>%
    transmute(
      imaging="resin LM",
      image_id=str_remove(name, "_c\\d+_[^_]+$"),
      uri, barcode=barcode_in(name), block_name=NA_character_, parse_name=name
    )
}

if ("tem_collection_table" %in% names(collection_tables)) {
  ## source_name: <grid>_<cut>_c001_<barcode|nobarcodeassigned>_<block name>[_P2]_correctionblend...
  image_tables$tem <- collection_tables$tem_collection_table %>%
    mutate(
      after_cell=str_remove(source_name, "^.*?_c\\d{3}_"),
      block_part=str_remove(after_cell, "^(\\d{9}|nobarcodeassigned)_") %>%
        str_remove("_correctionblend.*$") %>%
        str_remove("\\.zarr$")
    ) %>%
    transmute(
      imaging="TEM",
      image_id=source_name,
      uri, barcode=str_match(after_cell, "^(\\d{9})_")[, 2],
      block_name=str_remove(block_part, "_P\\d+$"),
      parse_name=block_part
    )
}

if ("hitt_collection_table" %in% names(collection_tables)) {
  image_tables$hitt <- collection_tables$hitt_collection_table %>%
    transmute(
      imaging="HiTT",
      image_id=name,
      uri, barcode=barcode_in(name),
      block_name=str_remove(name, "_P\\d+$"),
      parse_name=name
    )
}

if ("sem_collection_table" %in% names(collection_tables)) {
  ## no sample name: site, date and time of day from the columns
  image_tables$sem <- collection_tables$sem_collection_table %>%
    transmute(
      imaging="SEM",
      image_id=name,
      uri, barcode=NA_character_, block_name=NA_character_,
      ## Date / Time / Size fraction from the image log: complete, sem_date has typos (Porto 20230102)
      fraction=na_if(`Size fraction`, "NA"),
      parse_name=paste(site, Date, na_if(Time, "NA"), sep="_")
    )
}

## EXM: the collection table names are <site code>_Img<image number>; sample (barcode of the chem fixed
## FA2 aliquot), sampling date / time and fraction come from the ExM tab of the annotation sheet
exm_annotation_file <- file.path(exm_annotation_dir, "ExM.tsv")
if ("exm_collection_table" %in% names(collection_tables) && file.exists(exm_annotation_file)) {
  ## the first row of the tab holds the column groups, the column names are in the second
  exm_annotations <- read_tsv(exm_annotation_file, skip=1, col_types=cols(.default=col_character())) %>%
    transmute(
      exm_barcode=na_if(Barcode, "NA"),
      site_code=unname(set_names(names(site_codes), site_codes)[str_remove(Supersite, ",.*$")]),
      image_no=as.integer(str_remove_all(`Image number`, "\\D")),
      sampling_date=ymd(`Sampling Date`, quiet=TRUE),
      sampling_time=na_if(`Sampling Time`, "NA"),
      fraction=str_replace(na_if(Fraction, "NA"), "-", "to"),
      fraction_imaged=na_if(`Fraction Imaged`, "NA")
    ) %>%
    filter(!is.na(site_code), !is.na(image_no)) %>%
    mutate(image_id=sprintf("%s_Img%04d", site_code, image_no)) %>%
    distinct(image_id, .keep_all=TRUE)

  image_tables$exm <- collection_tables$exm_collection_table %>%
    one_per_uri() %>%
    transmute(imaging="ExM", image_id=str_remove(name, "_c\\d+_[^_]+$"), uri) %>%
    left_join(exm_annotations, by="image_id") %>%
    transmute(
      imaging, image_id, uri,
      barcode=NA_character_, block_name=NA_character_,
      exm_barcode,
      fraction,
      parse_name=paste(site_code, format(sampling_date, "%Y%m%d"), sampling_time, sep="_")
    )
  message("EXM images without annotation: ", sum(is.na(image_tables$exm$exm_barcode)))
}

imaging_map <- image_tables %>%
  bind_rows() %>%
  link_images()

## SEM: images of the chemically fixed samples, matched by site, date, time of day and size fraction
## (several tubes of one fraction: the first one), else by site, date and time of day only
chemfix_groups <- function(...) {
  chemfix_events %>%
    filter(!is.na(date), !is.na(event_id)) %>%
    group_by(site, date, time_of_day, ...) %>%
    summarise(cf_barcode=first(barcode), cf_event_id=first(event_id), .groups="drop")
}
sem_by_fraction <- imaging_map %>%
  filter(imaging == "SEM") %>%
  select(image_id, name_site, name_date, name_time, fraction) %>%
  inner_join(chemfix_groups(fraction), by=c(name_site="site", name_date="date", name_time="time_of_day", fraction="fraction")) %>%
  mutate(sem_link="sem_chemfix")
sem_by_slot <- imaging_map %>%
  filter(imaging == "SEM", !image_id %in% sem_by_fraction$image_id) %>%
  select(image_id, name_site, name_date, name_time) %>%
  inner_join(chemfix_groups(), by=c(name_site="site", name_date="date", name_time="time_of_day")) %>%
  mutate(sem_link="sem_chemfix_other_fraction")

imaging_map <- imaging_map %>%
  left_join(
    bind_rows(sem_by_fraction, sem_by_slot) %>% select(image_id, cf_barcode, cf_event_id, sem_link),
    by="image_id"
  ) %>%
  mutate(
    link=coalesce(sem_link, link),
    sample_barcode=coalesce(cf_barcode, sample_barcode),
    sample_layer=if_else(!is.na(cf_barcode), "Chem fixation", sample_layer),
    event_id=coalesce(cf_event_id, event_id),
    ## EXM: chemically fixed FA2 aliquots (own barcodes, not in the chem fixation sheet), event by site / date / time
    sample_barcode=if_else(imaging == "ExM", exm_barcode, sample_barcode),
    sample_layer=if_else(imaging == "ExM" & !is.na(event_id), "Chem fixation", sample_layer),
    link=if_else(imaging == "ExM" & !is.na(event_id), "exm_annotation", link),
    ## TEM / HiTT are always done on HPF resin blocks: blocks missing in the HPF sheet still count as HPF
    sample_layer=if_else(imaging %in% c("TEM", "HiTT") & is.na(sample_layer) & !is.na(event_id), "HPF", sample_layer),
    ## sample identity for counting samples: barcode, else the block name
    sample_key=coalesce(sample_barcode, block_name)
  ) %>%
  select(-cf_barcode, -cf_event_id, -sem_link, -exm_barcode)

imaging_map %>% count(imaging, link) %>% print(n=Inf)
imaging_map %>%
  filter(link == "none") %>%
  count(imaging, name_site, name_date, name_time) %>%
  print(n=Inf)

write_tsv(imaging_map, file=file.path(ov_dir, "imaging_map.tsv"))


## ---- events without any sample of a layer ----
event_coverage <- events %>%
  select(event_id, site, date, event_time, event_site_id) %>%
  mutate(
    n_hpf=map_int(event_id, ~sum(hpf_events$event_id %in% .x)),
    n_cryo_pf=map_int(event_id, ~sum(cryo_pf_events$event_id %in% .x)),
    n_cryo_wm=map_int(event_id, ~sum(cryo_wm_events$event_id %in% .x)),
    n_metab=map_int(event_id, ~sum(metab_events$event_id %in% .x)),
    n_chemfix=map_int(event_id, ~sum(chemfix_events$event_id %in% .x))
  ) %>%
  ## images per imaging modality and event
  left_join(
    imaging_map %>%
      filter(!is.na(event_id)) %>%
      count(event_id, imaging) %>%
      mutate(imaging=paste0("n_img_", str_replace_all(str_to_lower(imaging), "[^a-z]+", "_"))) %>%
      pivot_wider(names_from=imaging, values_from=n, values_fill=0),
    by="event_id"
  ) %>%
  mutate(across(starts_with("n_img_"), ~coalesce(.x, 0L)))
event_coverage %>%
  group_by(site) %>%
  summarise(events=n(), with_hpf=sum(n_hpf > 0), with_cryo_pf=sum(n_cryo_pf > 0), with_cryo_wm=sum(n_cryo_wm > 0), with_metab=sum(n_metab > 0), with_chemfix=sum(n_chemfix > 0))

write_tsv(event_coverage, file=file.path(ov_dir, "event_coverage.tsv"))
