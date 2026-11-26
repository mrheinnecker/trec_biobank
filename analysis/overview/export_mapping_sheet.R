## Biobank overview: all links between the layers of the overview figure, one tab per link, written
## to a Google sheet (and as tsv into <snapshot>/mapping_sheet/). Reads a data snapshot of
## biobank_overview.R, so the sheet shows exactly what the figures use.
##
## tabs (from the inside to the outside of the radial chain figure):
##   1_event_to_raw_sample     every raw sample (HPF, cryo PF, cryo WM, MetaB, Chem fixation) -> sampling event
##   2_raw_sample_to_derived   HPF sample -> resin block (Epon / Lowicryl) / kept in cryo / discarded;
##                             MetaB sample -> V4V5 rRNA sequencing (one per sample)
##   3_to_overview_scan        HiTT / resin LM / cryo LM datasets -> sample (barcode) or resin block
##   4_to_images               TEM -> HiTT scan of the same block or the block directly, SEM -> chem fixed
##                             tube, ExM -> chem fixed FA2 aliquot
##   summary                   rows per tab and status
## status: clear (unique match), guessed (several candidates, one picked by a rule: check), dropped
## (no match: not in the figure); in_figure = counted in the figure (status not dropped, site in the figure)
##
## env (set by analysis/launcher.sh): OVERVIEW_DATA_SNAPSHOT (folder name in <overview_root>/data, empty = newest),
## TREC_GOOGLE_KEY (service-account key with edit access to the sheet), OVERVIEW_MAPPING_SHEET (sheet URL,
## empty = only the tsv files)

library(tidyverse)

## data locations: analysis/paths.R (<data_root>/analysis/overview)
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
  dir <- if (is.na(f)) file.path(if (.Platform$OS.type == "windows") "C:/repos" else "/mnt/c/repos", "trec_biobank/analysis/overview") else dirname(normalizePath(f))
  source(file.path(dir, "..", "paths.R"), local=globalenv())
})
overview_root <- trec_paths$overview
mapping_sheet <- Sys.getenv("OVERVIEW_MAPPING_SHEET", "https://docs.google.com/spreadsheets/d/1DabzZ1275Z4Agh_ertKYsTSQXndRIsf-qnUz4VQvj2E/edit")
google_key <- Sys.getenv("TREC_GOOGLE_KEY", "")
data_snapshot <- Sys.getenv("OVERVIEW_DATA_SNAPSHOT", "")

ov_dir <- if (data_snapshot == "") {
  sort(list.dirs(file.path(overview_root, "data"), recursive=FALSE), decreasing=TRUE)[1]
} else {
  file.path(overview_root, "data", data_snapshot)
}
message("Data snapshot: ", ov_dir)
out_dir <- file.path(ov_dir, "mapping_sheet")
dir.create(out_dir, showWarnings=FALSE)

read_chr <- function(f) read_tsv(file.path(ov_dir, f), col_types=cols(.default=col_character()))

## sites shown in the figure (site_palette of biobank_overview_figure.R)
figure_sites <- c("Villefranche", "Iceland", "Roscoff", "Tallinn", "Kristineberg", "Bilbao", "Vigo", "Porto", "Barcelona", "Naples", "Athens")

events <- read_chr("event_coverage.tsv") %>% select(event_id, event_site=site, event_date=date, event_time, event_site_id)
overrides <- if (file.exists(file.path(ov_dir, "event_overrides.tsv"))) {
  read_tsv(file.path(ov_dir, "event_overrides.tsv"), col_types=cols(.default=col_character()), comment="#") %>%
    filter(!is.na(id)) %>%
    select(table, id, override_note=note)
} else {
  tibble(table=character(), id=character(), override_note=character())
}

sample_status <- function(q) case_when(
  q %in% c("date", "date_ampm", "event_id", "manual", "pilot_slot", "pilot_date") ~ "clear",
  q %in% c("untimed_event", "nearest_time") ~ "guessed",
  TRUE ~ "dropped"
)
how_sample <- c(
  date="only event at the site that day",
  date_ampm="only event in the same half day (AM / MID / PM)",
  untimed_event="several events, none in that half day: the one without collection time",
  nearest_time="several events: the closest to AM 9:00 / MID 12:00 / PM 15:00",
  no_sample_time="several events that day, sample has no time of day",
  no_event="no event at the site on that date",
  event_id="event ID in the sample log (Vigo log)",
  manual="set by hand (event_overrides.tsv)",
  manual_none="taken out by hand (event_overrides.tsv)",
  pilot_slot="pilots log: event with the same date and time of day",
  pilot_date="pilots log: sample without time of day, only event at the site that day",
  pilot_unclear="pilots log: no clear event (time of day missing or not in the log), dropped (to check)"
)


## ---- sort order ----
## certainty: 1 = least certain (no link) ... higher = more certain; rows are sorted
##   1. dropped (no link), 2. linked but not in the figure (e.g. site not shown, sample without event),
##   3. in the figure; within these by certainty (least certain first), then by type and site
certainty_sample <- c(no_event=1, no_sample_time=1, manual_none=1, nearest_time=2, pilot_unclear=1, untimed_event=3,
                      date_ampm=4, date=5, pilot_date=5, event_id=6, pilot_slot=6, manual=7)
certainty_image <- c(none=1, manual_none=1, name_event=2, sem_chemfix_other_fraction=3, sem_chemfix=4,
                     exm_annotation=5, block_name=6, barcode=7, manual=8)
sort_tab <- function(tab, ...) {
  tab %>%
    arrange(case_when(status == "dropped" ~ 1, !in_figure ~ 2, TRUE ~ 3), certainty, ...) %>%
    relocate(certainty, .after=status)
}


## ---- 1: sampling event -> raw sample ----
sample_files <- c(
  HPF="hpf_event_map.tsv", `cryo PF`="cryo_pf_event_map.tsv", `cryo WM`="cryo_wm_event_map.tsv",
  MetaB="metab_event_map.tsv", `Chem fixation`="chemfix_event_map.tsv"
)
samples <- imap_dfr(sample_files, ~mutate(read_chr(.x), sample_type=.y))

tab1 <- samples %>%
  left_join(events %>% select(event_id, event_site), by="event_id") %>%
  left_join(overrides, by=c("sample_type"="table", "barcode"="id")) %>%
  transmute(
    sample_type, barcode,
    sheet_site=site, sheet_date=date_raw, sheet_time_of_day=time_of_day,
    status=sample_status(match_quality),
    certainty=unname(certainty_sample[match_quality]),
    in_figure=status != "dropped" & !is.na(event_id) & event_site %in% figure_sites,
    match_quality, how_linked=unname(how_sample[match_quality]),
    event_site_id, event_time, event_id,
    candidate_events, candidate_event_ids, override_note
  ) %>%
  sort_tab(sample_type, sheet_site, sheet_date)


## ---- 2: raw sample -> derived sample ----
hpf <- read_chr("hpf_event_map.tsv")
images <- read_chr("imaging_map.tsv") %>%
  left_join(events %>% select(event_id, event_site), by="event_id")
resin_of <- function(x) if_else(str_detect(x, regex("_low", ignore_case=TRUE)), "Lowicryl", "Epon")

## resin blocks linked to the event only (Vigo log)
blocks_extra <- if (file.exists(file.path(ov_dir, "resin_blocks_extra.tsv"))) read_chr("resin_blocks_extra.tsv") else
  tibble(site=character(), event_id=character(), cryo_hpf_id=character(), resin=character(), block_name=character(), sample_key=character(),
         how_linked=character(), match_quality=character())

## in_figure of the sample from its own event (not by barcode: some HPF rows have no barcode);
## the Vigo carriers are not linked to single blocks: their blocks are listed per HPF puck position below
hpf_derived <- hpf %>%
  filter(site != "Vigo" | is.na(site)) %>%
  left_join(events %>% select(event_id, event_site), by="event_id") %>%
  mutate(sample_in_figure=!is.na(event_id) & event_site %in% figure_sites) %>%
  transmute(
    from_type="HPF", from_id=barcode, site,
    derived=case_when(
      block_name == "discarded" ~ "discarded",
      block_name == "block?" ~ "resin block (name unclear)",
      !is.na(block_name) ~ resin_of(block_name),
      !is.na(cryo_name) ~ "kept in cryo",
      TRUE ~ "nothing recorded"
    ),
    derived_id=coalesce(if_else(block_name %in% c("discarded"), NA_character_, block_name), cryo_name),
    status=case_when(
      derived == "nothing recorded" ~ "dropped",
      derived == "resin block (name unclear)" ~ "guessed",
      TRUE ~ "clear"
    ),
    how_linked=case_when(
      derived %in% c("Epon", "Lowicryl") ~ "block name in the HPF sheet (_low = Lowicryl, else Epon)",
      derived == "resin block (name unclear)" ~ "block name 'block?' in the HPF sheet: counted as Epon",
      derived == "kept in cryo" ~ "cryo name in the HPF sheet",
      derived == "discarded" ~ "block name 'discarded' in the HPF sheet: not in the figure",
      TRUE ~ "no block or cryo name in the HPF sheet"
    ),
    in_figure=coalesce(sample_in_figure, FALSE) & status != "dropped" & derived != "discarded",
    certainty=case_when(status == "dropped" ~ 1, derived == "resin block (name unclear)" ~ 3, derived == "discarded" ~ 5, TRUE ~ 6)
  )

## blocks only known from TEM / HiTT / resin LM images (no block name in the HPF sheet)
image_blocks <- images %>%
  filter(imaging %in% c("TEM", "HiTT", "resin LM"), sample_layer == "HPF") %>%
  mutate(block=coalesce(block_name, image_id)) %>%
  filter(!sample_key %in% (hpf %>% filter(!is.na(block_name), block_name != "discarded") %>% pull(barcode))) %>%
  filter(!sample_key %in% blocks_extra$sample_key) %>%
  group_by(sample_key) %>%
  summarise(
    site=first(event_site), block=first(block), n_images=n(),
    imaging=paste(sort(unique(imaging)), collapse=" + "), link=paste(sort(unique(link)), collapse=" + "),
    .groups="drop"
  ) %>%
  transmute(
    from_type="HPF", from_id=sample_key, site,
    derived=resin_of(block), derived_id=block,
    status=if_else(str_detect(link, "name_event"), "guessed", "clear"),
    how_linked=paste0("block not in the HPF sheet, known from ", n_images, " ", imaging, " image(s) (linked by ", link, ")"),
    in_figure=site %in% figure_sites,
    certainty=if_else(status == "guessed", 2, 4)
  )

metab_derived <- tab1 %>%
  filter(sample_type == "MetaB") %>%
  transmute(
    from_type="MetaB", from_id=barcode, site=sheet_site,
    derived="V4V5 rRNA sequencing", derived_id=barcode,
    status=if_else(in_figure, "clear", "dropped"),
    how_linked="one sequencing library per MetaB sample (assumed, not from a sheet)",
    in_figure,
    certainty=if_else(status == "dropped", 1, 6)
  )

event_blocks <- blocks_extra %>%
  left_join(events %>% select(event_id, event_site), by="event_id") %>%
  transmute(
    from_type="HPF", from_id=coalesce(cryo_hpf_id, event_id), site=coalesce(event_site, site),
    derived=resin, derived_id=sample_key,
    status=if_else(is.na(event_id), "dropped", "clear"),
    how_linked=paste0(how_linked, " (puck position with several carriers, not a single carrier)"),
    in_figure=!is.na(event_id) & site %in% figure_sites,
    certainty=if_else(is.na(event_id), 1, 6)
  )

tab2 <- bind_rows(hpf_derived, image_blocks, event_blocks, metab_derived) %>%
  sort_tab(from_type, site, from_id)


## ---- 3 / 4: images ----
image_status <- function(link) case_when(
  link %in% c("barcode", "block_name", "exm_annotation", "sem_chemfix", "manual") ~ "clear",
  link %in% c("name_event", "sem_chemfix_other_fraction") ~ "guessed",
  TRUE ~ "dropped"
)
how_image <- c(
  barcode="barcode in the image name -> sample",
  block_name="block name in the image name -> HPF sample with that block",
  name_event="site / date / time of day in the name -> sampling event only, sample unknown",
  none="nothing found in the name",
  sem_chemfix="site / date / time / size fraction -> first chem fixed tube of that slot",
  sem_chemfix_other_fraction="site / date / time (no tube of that fraction) -> first chem fixed tube of that slot",
  exm_annotation="ExM annotation tab: FA2 aliquot barcode, event by site / date / time",
  manual="set by hand (event_overrides.tsv)",
  manual_none="taken out by hand (event_overrides.tsv)"
)
image_rows <- images %>%
  left_join(overrides %>% filter(table == "image") %>% select(-table), by=c("image_id"="id")) %>%
  mutate(
    status=image_status(link),
    certainty=unname(certainty_image[link]),
    in_figure=status != "dropped" & !is.na(event_id) & !is.na(sample_layer) & event_site %in% figure_sites,
    how_linked=unname(how_image[link])
  )

tab3 <- image_rows %>%
  filter(imaging %in% c("HiTT", "resin LM", "cryo LM PF", "cryo LM WM", "cryo LM HPF")) %>%
  transmute(
    imaging, image_id, status, certainty, in_figure, link, how_linked,
    sample_layer, sample_barcode, block_name, resin=if_else(imaging %in% c("HiTT", "resin LM"), resin_of(coalesce(block_name, image_id)), NA_character_),
    site=event_site, event_id, name_site, name_date, name_time, override_note, uri
  ) %>%
  sort_tab(imaging, site, image_id)

## TEM comes from HiTT when the same sample / block has HiTT images
hitt_keys <- images %>% filter(imaging == "HiTT", !is.na(sample_key)) %>% pull(sample_key) %>% unique()
tab4 <- image_rows %>%
  filter(imaging %in% c("TEM", "SEM", "ExM")) %>%
  transmute(
    imaging, image_id, status, certainty, in_figure, link, how_linked,
    parent=case_when(
      imaging == "TEM" & sample_key %in% hitt_keys ~ "HiTT scan of the same block",
      imaging == "TEM" ~ paste0("resin block directly (", resin_of(coalesce(block_name, image_id)), ")"),
      imaging == "SEM" ~ "chem fixed tube",
      imaging == "ExM" ~ "chem fixed FA2 aliquot"
    ),
    parent_id=coalesce(block_name, sample_barcode),
    sample_layer, sample_barcode, block_name, sample_key,
    site=event_site, event_id, name_site, name_date, name_time, fraction, override_note, uri
  ) %>%
  sort_tab(imaging, site, image_id)


## ---- summary + write ----
tabs <- list(
  `1_event_to_raw_sample`=tab1,
  `2_raw_sample_to_derived`=tab2,
  `3_to_overview_scan`=tab3,
  `4_to_images`=tab4
)
type_col <- c(`1_event_to_raw_sample`="sample_type", `2_raw_sample_to_derived`="from_type",
              `3_to_overview_scan`="imaging", `4_to_images`="imaging")
summary_tab <- imap_dfr(tabs, function(tab, name) {
  tab %>%
    count(type=.data[[type_col[[name]]]], status, in_figure) %>%
    mutate(tab=name, .before=1)
}) %>%
  pivot_wider(names_from=c(status, in_figure), values_from=n, values_fill=0, names_glue="{status}{if_else(in_figure, '', ' (not in figure)')}") %>%
  mutate(total=rowSums(across(where(is.numeric)))) %>%
  mutate(snapshot=basename(ov_dir), written=format(Sys.time(), "%Y-%m-%d %H:%M"))
summary_tab <- summary_tab %>% select(tab, type, any_of(c("clear", "guessed", "dropped")), everything())
tabs$summary <- summary_tab
print(summary_tab, n=Inf, width=Inf)

iwalk(tabs, ~write_tsv(.x, file.path(out_dir, paste0(.y, ".tsv")), na=""))
message("tsv files: ", out_dir)

if (mapping_sheet != "") {
  if (google_key != "") googlesheets4::gs4_auth(path=google_key)
  iwalk(tabs, function(tab, name) {
    googlesheets4::sheet_write(mutate(tab, across(everything(), as.character)), ss=mapping_sheet, sheet=name)
    message("sheet tab ", name, ": ", nrow(tab), " rows")
  })
}
