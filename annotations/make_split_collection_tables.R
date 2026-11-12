## Distribute TEM images to the annotators and collect their annotations.
##
## Every image of the TEM collection table is annotated twice by different people:
##   round 1 (annotation): name in annotated_by
##   round 2 (validation): name in validated_by
##
## Sheets in the annotation spreadsheet (--collection_table_url):
##   <annotations_sheet> (main): one row per finished round (round, annotated_by = who did it)
##   <status_sheet>:             one row per image: annotated_by (round 1), validated_by (round 2)
##   <person> tabs:              the current assignment of each annotator
##
## Every run:
##   1. reads main, status and all personal tabs and backs them up to <annotation_log_dir>/<timestamp>/
##   2. moves finished rows (values + name filled) from the personal tabs into main
##   3. keeps started but unsigned rows (values, no name) in the same personal tab
##   4. rewrites status and refills the personal tabs: all open images are split equally
##      over the assignees, first annotations before validations
##      (--max_rows_per_person optionally caps the number of rows per person)
##
## People who leave: move them from --assignees to --former_assignees. Their tab is still
## read and backed up, finished rows go to main, then the tab is emptied.
##
## The main table must already be in the round format (see annotations/migrate_annotation_rounds.R).

library(tidyverse)
library(getopt)

spec <- matrix(c(
  "collection_table_url", "u", 1, "character",
  "source_collection_table_url", "c", 1, "character",
  "source_collection_sheet", "t", 1, "character",
  "local_source_collection_table", "f", 1, "character",
  "google_key", "k", 1, "character",
  "sheet_mode", "m", 1, "character",
  "local_outdir", "o", 1, "character",
  "annotation_log_dir", "a", 1, "character",
  "annotations_sheet", "n", 1, "character",
  "status_sheet", "s", 1, "character",
  "assignees", "g", 1, "character",
  "former_assignees", "r", 1, "character",
  "max_rows_per_person", "w", 1, "character"
), ncol=4, byrow=TRUE)
opt <- getopt(spec)

# Debug helper: uncomment to run interactively without command-line parsing.
# opt <- list(
#   collection_table_url = "https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=0#gid=0",
#   source_collection_table_url = "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=2014112369#gid=2014112369",
#   source_collection_sheet = "tem_collection_table",
#   google_key = "/g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json",
#   sheet_mode = "google",
#   annotation_log_dir = "/g/schwab/tem_screen/annotations/log",
#   annotations_sheet = "main",
#   status_sheet = "status",
#   assignees = "marco,chandni,yannick,karel,ella",
#   former_assignees = "viktoria",
#   max_rows_per_person = NA
# )

arg_file <- commandArgs(FALSE) %>%
  .[str_detect(., "^--file=")] %>%
  str_remove("^--file=") %>%
  .[1]

script_dir <- if (!is.na(arg_file)) {
  dirname(normalizePath(arg_file))
} else {
  getwd()
}

source(file.path(script_dir, "annotation_rounds_common.R"))

opt_or <- function(value, default) {
  if (is.null(value) || is.na(value)) default else value
}

collection_table_url <- opt_or(
  opt$collection_table_url,
  "https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=0#gid=0"
)
source_collection_table_url <- opt_or(
  opt$source_collection_table_url,
  "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=2014112369#gid=2014112369"
)
source_collection_sheet <- opt_or(opt$source_collection_sheet, "tem_collection_table")
local_source_collection_table <- opt_or(opt$local_source_collection_table, "tem_collection_table.tsv")
google_key <- opt_or(opt$google_key, file.path(script_dir, "..", "trec-tem-screen-e98a2e03f58b.json"))
sheet_mode <- opt_or(opt$sheet_mode, "google")
local_outdir <- opt_or(opt$local_outdir, "split_collection_tables")
annotation_log_dir <- opt_or(opt$annotation_log_dir, "/g/schwab/tem_screen/annotations/log")
annotations_sheet <- opt_or(opt$annotations_sheet, "main")
status_sheet <- opt_or(opt$status_sheet, "status")

parse_names <- function(x) {
  str_split(x, ",")[[1]] %>%
    str_trim() %>%
    str_to_lower() %>%
    discard(~.x == "" | is.na(.x)) %>%
    unique()
}
people <- parse_names(opt_or(opt$assignees, "marco,chandni,yannick,karel,ella"))
former_people <- parse_names(opt_or(opt$former_assignees, ""))
if (length(people) == 0) {
  stop("--assignees must contain at least one sheet/person name.")
}
if (length(intersect(people, former_people)) > 0) {
  stop("Names in both --assignees and --former_assignees: ",
       paste(intersect(people, former_people), collapse=", "))
}
if (any(c(annotations_sheet, status_sheet) %in% c(people, former_people))) {
  stop("An assignee name is the same as the main or status sheet name.")
}
all_tabs <- c(people, former_people)

max_rows_per_person <- opt$max_rows_per_person
if (is.null(max_rows_per_person) || is.na(max_rows_per_person)) {
  max_rows_per_person <- rep(Inf, length(people))
} else {
  max_rows_per_person <- str_split(max_rows_per_person, ",", simplify=FALSE)[[1]] %>%
    str_trim() %>%
    as.numeric()
  if (any(is.na(max_rows_per_person))) {
    stop("--max_rows_per_person must be one number or a comma-separated number for each assignee.")
  }
  if (length(max_rows_per_person) == 1) {
    max_rows_per_person <- rep(max_rows_per_person, length(people))
  } else if (length(max_rows_per_person) != length(people)) {
    stop(
      "--max_rows_per_person must contain either one value or exactly ",
      length(people), " values (one per assignee)."
    )
  }
}
if (any(max_rows_per_person < 0)) {
  stop("--max_rows_per_person values must be 0 or greater.")
}
max_rows_per_person <- set_names(max_rows_per_person, people)

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

## ---- google / local io ----

google_auth <- function() {
  library(googlesheets4)
  library(googledrive)
  if (is.null(google_key) || is.na(google_key) || !file.exists(google_key)) {
    stop("--google_key is required and must exist when --sheet_mode google")
  }
  gs4_auth(path=google_key)
  drive_auth(path=google_key)
}

freeze_first_column <- function(ss, sheet_name) {
  spreadsheet_id <- as.character(googlesheets4::as_sheets_id(ss))
  sheet_id <- googlesheets4::sheet_properties(ss) %>%
    filter(name == sheet_name) %>%
    pull(id)

  if (length(sheet_id) != 1) {
    warning("Could not freeze first column for sheet: ", sheet_name)
    return(invisible(NULL))
  }

  response <- httr::POST(
    url=sprintf(
      "https://sheets.googleapis.com/v4/spreadsheets/%s:batchUpdate",
      spreadsheet_id
    ),
    config=googlesheets4::gs4_token(),
    body=list(
      requests=list(
        list(
          updateSheetProperties=list(
            properties=list(
              sheetId=sheet_id,
              gridProperties=list(frozenColumnCount=1)
            ),
            fields="gridProperties.frozenColumnCount"
          )
        )
      )
    ),
    encode="json"
  )

  httr::stop_for_status(response)
  invisible(NULL)
}

sheet_names <- NULL
read_google_sheet_if_exists <- function(ss, sheet_name) {
  if (is.null(sheet_names)) {
    sheet_names <<- googlesheets4::sheet_properties(ss)$name
  }
  if (!sheet_name %in% sheet_names) {
    warning("Sheet does not exist yet, treating as empty: ", sheet_name)
    return(tibble())
  }
  googlesheets4::read_sheet(ss, sheet=sheet_name, col_types="c")
}

read_local_if_exists <- function(sheet_name) {
  file <- file.path(local_outdir, paste0(sheet_name, ".tsv"))
  if (file.exists(file)) {
    read_tsv(file, col_types=cols(.default=col_character()))
  } else {
    tibble()
  }
}

read_annotation_sheet <- function(sheet_name) {
  if (sheet_mode == "google") {
    read_google_sheet_if_exists(collection_table_url, sheet_name)
  } else {
    read_local_if_exists(sheet_name)
  }
}

write_annotation_sheet <- function(table, sheet_name, freeze=FALSE) {
  if (sheet_mode == "google") {
    googlesheets4::write_sheet(table, ss=collection_table_url, sheet=sheet_name)
    if (freeze) {
      freeze_first_column(collection_table_url, sheet_name)
    }
  } else {
    dir.create(local_outdir, recursive=TRUE, showWarnings=FALSE)
    write_tsv(table, file.path(local_outdir, paste0(sheet_name, ".tsv")), na="")
  }
}

read_source_collection_table <- function() {
  if (sheet_mode == "google") {
    googlesheets4::read_sheet(source_collection_table_url, sheet=source_collection_sheet, col_types="c")
  } else if (file.exists(local_source_collection_table)) {
    read_tsv(local_source_collection_table, col_types=cols(.default=col_character()))
  } else {
    stop("Local source collection table does not exist: ", local_source_collection_table)
  }
}

standardize_collection_table <- function(col_table) {
  if (!"name" %in% names(col_table)) {
    stop("The source collection table needs a name column.")
  }
  if (!"source_name" %in% names(col_table)) {
    if ("uri" %in% names(col_table)) {
      col_table$source_name <- basename(str_remove(as.character(col_table$uri), "/$"))
    } else {
      col_table$source_name <- col_table$name
    }
  }
  if (!"site" %in% names(col_table)) {
    col_table$site <- str_extract(col_table$source_name, "ATH|BAR|KRI|TAL|NAP|BIL|POR")
  }
  if (!"cell_id" %in% names(col_table)) {
    col_table$cell_id <- str_extract(col_table$source_name, "c0\\d+")
  }
  if (!"size_frac" %in% names(col_table)) {
    col_table$size_frac <- str_extract(col_table$source_name, "\\d+to\\d+")
  }
  if (!"sampling_time" %in% names(col_table)) {
    col_table$sampling_time <- str_extract(col_table$source_name, "_(AM|PM|MID|TARA)_") %>%
      str_remove_all("_")
  }
  if (!"exclusive" %in% names(col_table)) {
    col_table$exclusive <- "TRUE"
  }

  col_table %>%
    filter(!is_blank(name)) %>%
    distinct(name, .keep_all=TRUE) %>%
    select(any_of(collection_identity_columns)) %>%
    arrange(site, name)
}

## ---- backup ----

write_backup_table <- function(table, path) {
  if (ncol(table) == 0) {
    write_lines(character(), path)
  } else {
    write_tsv(table, path, na="")
  }
}

## raw copy of every sheet before anything is changed; stops if incomplete
backup_annotations <- function(sources) {
  backup_dir <- file.path(annotation_log_dir, run_timestamp)
  if (dir.exists(backup_dir)) {
    stop("Backup directory already exists; refusing to overwrite it: ", backup_dir)
  }
  dir.create(backup_dir, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(backup_dir)) {
    stop("Could not create backup directory: ", backup_dir)
  }

  write_backup_table(sources$main, file.path(backup_dir, paste0(annotations_sheet, ".tsv")))
  write_backup_table(sources$status, file.path(backup_dir, paste0(status_sheet, ".tsv")))
  iwalk(sources$assignments, function(table, person) {
    write_backup_table(table, file.path(backup_dir, paste0(person, ".tsv")))
  })
  write_tsv(
    tibble(
      setting=c(
        "sheet_mode", "collection_table_url", "annotations_sheet", "status_sheet",
        "source_collection_table_url", "source_collection_sheet",
        "assignees", "former_assignees", "max_rows_per_person"
      ),
      value=c(
        sheet_mode, collection_table_url, annotations_sheet, status_sheet,
        source_collection_table_url, source_collection_sheet,
        paste(people, collapse=","), paste(former_people, collapse=","),
        paste(max_rows_per_person, collapse=",")
      )
    ),
    file.path(backup_dir, "run_configuration.tsv")
  )

  expected_files <- c(
    paste0(c(annotations_sheet, status_sheet, all_tabs), ".tsv"),
    "run_configuration.tsv"
  )
  missing_files <- expected_files[!file.exists(file.path(backup_dir, expected_files))]
  if (length(missing_files) > 0) {
    stop("Backup is incomplete; no sheet will be overwritten. Missing: ",
         paste(missing_files, collapse=", "))
  }

  message("Complete pre-write backup saved to: ", backup_dir)
  backup_dir
}

## ---- reading the current state ----

read_current_annotation_sources <- function() {
  list(
    main=read_annotation_sheet(annotations_sheet),
    status=read_annotation_sheet(status_sheet),
    assignments=set_names(map(all_tabs, read_annotation_sheet), all_tabs)
  )
}

parse_main_table <- function(main) {
  if (nrow(main) == 0) {
    return(tibble(name=character(), round=integer(), annotated_by=character()))
  }
  if (!"round" %in% names(main) || !"annotated_by" %in% names(main)) {
    stop(
      "The '", annotations_sheet, "' sheet has no round/annotated_by columns (old format). ",
      "Run annotations/migrate_annotation_rounds.R first."
    )
  }
  main <- main %>%
    filter(!is_blank(name)) %>%
    mutate(
      round=suppressWarnings(as.integer(round)),
      annotated_by=normalize_person(annotated_by)
    )
  bad <- main %>% filter(!round %in% c(1L, 2L) | is.na(annotated_by))
  if (nrow(bad) > 0) {
    stop("Rows in '", annotations_sheet, "' without a valid round (1/2) or annotated_by: ",
         paste(head(bad$name, 20), collapse=", "))
  }
  dup <- main %>% count(name, round) %>% filter(n > 1)
  if (nrow(dup) > 0) {
    stop("Images with the same round twice in '", annotations_sheet, "': ",
         paste(head(dup$name, 20), collapse=", "))
  }
  same_person <- main %>% count(name, annotated_by) %>% filter(n > 1)
  if (nrow(same_person) > 0) {
    stop("Images annotated twice by the same person in '", annotations_sheet, "': ",
         paste(head(same_person$name, 20), collapse=", "))
  }
  main
}

collect_tab_rows <- function(assignments) {
  imap_dfr(assignments, function(table, person) {
    if (nrow(table) == 0 || !"name" %in% names(table)) {
      return(tibble())
    }
    table %>%
      mutate(across(everything(), as.character)) %>%
      mutate(.sheet=person, .run=run_timestamp, .source=paste0("tab ", person))
  })
}

## ---- assignment ----

## first annotations: contiguous blocks (sorted by site, name) per person, up to capacity
split_contiguous <- function(rows, capacity) {
  split_tables <- set_names(replicate(length(people), rows[0, ], simplify=FALSE), people)
  if (nrow(rows) == 0 || all(capacity == 0)) {
    return(split_tables)
  }

  rows <- rows %>% arrange(site, name)
  if (!any(is.infinite(capacity))) {
    rows <- rows %>% head(sum(capacity))
  }

  next_row <- 1L
  for (person_index in seq_along(people)) {
    remaining_rows <- nrow(rows) - next_row + 1L
    if (remaining_rows <= 0) {
      break
    }
    later <- if (person_index == length(people)) numeric() else capacity[(person_index + 1L):length(people)]
    later_capacity <- sum(later)
    ## only people who can still take rows share the remaining ones
    remaining_people <- sum(capacity[person_index] > 0) + sum(later > 0)
    if (remaining_people == 0) {
      break
    }
    minimum_needed_here <- if (is.infinite(later_capacity)) 0 else max(0, remaining_rows - later_capacity)
    balanced_share <- ceiling(remaining_rows / remaining_people)
    take <- min(capacity[[person_index]], max(minimum_needed_here, balanced_share), remaining_rows)
    if (take > 0) {
      last_row <- next_row + take - 1L
      split_tables[[people[[person_index]]]] <- rows[next_row:last_row, , drop=FALSE]
      next_row <- last_row + 1L
    }
  }

  split_tables
}

## validations: one by one to the person with the fewest rows who did not do round 1
assign_validations <- function(rows, split_tables, capacity, offset) {
  if (nrow(rows) == 0) {
    return(split_tables)
  }
  ## rows already in the tab (kept unsigned rows) count for the balance
  counts <- map_int(split_tables, nrow) + offset[people]
  rows <- rows %>% arrange(site, name)
  picked <- set_names(vector("list", length(people)), people)

  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop=FALSE]
    allowed <- setdiff(people, row$annotated_by)
    allowed <- allowed[counts[allowed] < capacity[allowed] + offset[allowed]]
    if (length(allowed) == 0) {
      next
    }
    target <- allowed[which.min(counts[allowed])]
    picked[[target]] <- c(picked[[target]], i)
    counts[[target]] <- counts[[target]] + 1L
  }

  imap(split_tables, function(table, person) {
    bind_rows(table, rows[picked[[person]], , drop=FALSE])
  })
}

## ---- main ----

if (sheet_mode == "google") {
  google_auth()
}

## 1. read everything and back it up before any change
annotation_sources <- read_current_annotation_sources()
backup_dir <- backup_annotations(annotation_sources)

collection <- read_source_collection_table() %>%
  standardize_collection_table()
if (nrow(collection) == 0) {
  stop("No rows found in source collection table: ", source_collection_sheet)
}

existing_main <- parse_main_table(annotation_sources$main)

## 2. classify the rows of the personal tabs
tab_rows <- classify_tab_rows(collect_tab_rows(annotation_sources$assignments))
if (nrow(tab_rows) == 0) {
  tab_rows <- tibble(name=character(), .status=character(), .person=character(), .sheet=character())
}
signed <- tab_rows %>% filter(.status == "signed")
in_progress <- tab_rows %>% filter(.status == "in_progress")

unknown <- signed %>% filter(!.person %in% all_tabs)
if (nrow(unknown) > 0) {
  stop(
    "Unknown annotator names in the personal tabs (fix the sheet or add them to name_aliases in ",
    "annotations/annotation_rounds_common.R). Nothing was written.\n",
    paste0("  tab ", unknown$.sheet, ": ", unknown$name, " -> '", coalesce(unknown$validated_by, unknown$annotated_by), "'",
           collapse="\n")
  )
}

## 3. finished rows -> main
resolved <- resolve_rounds(existing_main, signed)
new_main <- resolved$rounds %>%
  attach_identity(collection) %>%
  sort_rounds()
issues <- resolved$issues

not_in_collection <- setdiff(unique(new_main$name), collection$name)
if (length(not_in_collection) > 0) {
  issues <- bind_rows(issues, tibble(
    type="not_in_collection_table_kept_in_main", name=not_in_collection, person=NA_character_,
    detail="image not in the TEM collection table"
  ))
}

done_round <- function(r) new_main %>% filter(round == r) %>% select(name, annotated_by)
round1 <- done_round(1L)
round2 <- done_round(2L)

## 4. started but unsigned rows stay with the tab they are in (if the round is still open);
##    in tabs of former assignees they go back to the pool
former_unsigned <- in_progress %>% filter(.sheet %in% former_people)
if (nrow(former_unsigned) > 0) {
  issues <- bind_rows(issues, former_unsigned %>% transmute(
    type="unsigned_values_of_former_assignee_returned_to_pool", name, person=.sheet,
    detail=paste0(.n_filled, " fields filled but not signed; values only in the backup")
  ))
}
in_progress <- in_progress %>% filter(.sheet %in% people)
if (nrow(in_progress) > 0) {
  in_progress <- in_progress %>%
    mutate(
      round=case_when(
        .tab_round %in% c(1L, 2L) ~ .tab_round,
        name %in% round1$name ~ 2L,
        TRUE ~ 1L
      ),
      .open=case_when(
        round == 1L ~ !name %in% round1$name,
        TRUE ~ !name %in% round2$name
      )
    ) %>%
    distinct(name, .keep_all=TRUE)
  closed <- in_progress %>% filter(!.open)
  if (nrow(closed) > 0) {
    issues <- bind_rows(issues, closed %>% transmute(
      type="unsigned_values_round_already_done_removed_from_tab", name, person=.sheet,
      detail=paste0("round ", round, " was finished by someone else; values only in the backup")
    ))
  }
  in_progress <- in_progress %>% filter(.open)
  if (nrow(in_progress) > 0) {
    issues <- bind_rows(issues, in_progress %>% transmute(
      type="unsigned_values_kept_in_tab", name, person=.sheet,
      detail=paste0(.n_filled, " fields filled but no name in ",
                    if_else(round == 1L, "annotated_by", "validated_by"))
    ))
  }
}

kept <- set_names(map(people, function(person) {
  rows <- in_progress %>% filter(.sheet == person)
  if (nrow(rows) == 0) {
    return(tibble())
  }
  rows %>%
    select(-any_of(c("annotated_by", "validated_by"))) %>%
    left_join(round1 %>% rename(.r1=annotated_by), by="name") %>%
    mutate(
      annotated_by=if_else(round == 2L, .r1, NA_character_),
      validated_by=NA_character_
    ) %>%
    select(-starts_with(".")) %>%
    select(-any_of(setdiff(collection_identity_columns, "name"))) %>%
    left_join(collection, by="name")
}), people)

busy <- unlist(map(kept, ~if ("name" %in% names(.x)) .x$name else character()))
capacity <- set_names(pmax(0, max_rows_per_person - map_int(kept, nrow)), people)

## 5. new assignments
first_pool <- collection %>%
  filter(!name %in% round1$name, !name %in% busy) %>%
  mutate(round=1L, annotated_by=NA_character_, validated_by=NA_character_)

second_pool <- collection %>%
  filter(name %in% round1$name, !name %in% round2$name, !name %in% busy) %>%
  left_join(round1, by="name") %>%
  mutate(round=2L, validated_by=NA_character_)

split_tables <- split_contiguous(first_pool, capacity)
split_tables <- assign_validations(second_pool, split_tables, capacity, map_int(kept, nrow))

blank_values <- set_names(rep(list(NA_character_), nrow(annotation_columns)), annotation_columns$short)
split_tables <- imap(split_tables, function(table, person) {
  new_rows <- table %>%
    mutate(round=as.integer(round)) %>%
    bind_cols(as_tibble(blank_values)[rep(1, nrow(table)), ])
  kept_rows <- kept[[person]]
  if (nrow(kept_rows) > 0) {
    kept_rows <- kept_rows %>% mutate(round=as.integer(round))
  }
  bind_rows(kept_rows, new_rows) %>%
    mutate(view=person, grid=person) %>%
    arrange(round, site, name) %>%
    order_tab_columns()
})
## tabs of former assignees are emptied (their content is in the backup)
former_tables <- set_names(map(former_people, ~split_tables[[1]][0, ]), former_people)

status <- make_status_table(collection, new_main)

## 6. save the new state into the backup folder, then write the sheets
write_backup_table(new_main, file.path(backup_dir, paste0("new_", annotations_sheet, ".tsv")))
write_backup_table(status, file.path(backup_dir, paste0("new_", status_sheet, ".tsv")))
iwalk(c(split_tables, former_tables), ~write_backup_table(.x, file.path(backup_dir, paste0("new_", .y, ".tsv"))))
if (nrow(issues) > 0) {
  write_tsv(issues, file.path(backup_dir, "issues.tsv"), na="")
}

write_annotation_sheet(new_main, annotations_sheet, freeze=TRUE)
write_annotation_sheet(status, status_sheet, freeze=TRUE)
iwalk(c(split_tables, former_tables), ~write_annotation_sheet(.x, .y, freeze=TRUE))

## ---- summary ----
n_added <- if (is.null(resolved$added)) 0L else nrow(resolved$added)
message("Rounds moved from the personal tabs into ", annotations_sheet, ": ", n_added)
message(
  "Images: ", nrow(collection), " | round 1 done: ", sum(collection$name %in% round1$name),
  " | round 2 done: ", sum(collection$name %in% round2$name)
)
if (nrow(issues) > 0) {
  message("Issues (also in ", file.path(backup_dir, "issues.tsv"), "):")
  print(issues %>% count(type), n=Inf)
}
print(tibble(
  sheet=names(split_tables),
  annotation=map_int(split_tables, ~sum(.x$round == 1L)),
  validation=map_int(split_tables, ~sum(.x$round == 2L)),
  unsigned_kept=map_int(kept, nrow)
))
