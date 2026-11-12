## Shared helpers for the TEM annotation rounds.
## Sourced by annotations/make_split_collection_tables.R and annotations/migrate_annotation_rounds.R.
##
## Every image is annotated twice by different people:
##   round 1 (annotation): the person writes their name into annotated_by
##   round 2 (validation): the person writes their name into validated_by
## The main table holds one row per finished round (round + annotated_by = who did it).

annotation_columns <- tribble(
  ~full, ~short,
  "life_status", "life_status",
  "major_group", "major_group",
  "taxomomic_class", "taxo_class",
  "cell_covering", "cell_cover",
  "thecal_struc", "thecal_struc",
  "nucleus",  "nucleus",
  "nucleolus",  "nucleolus",
  "golgi","golgi",
  "er", "ER",
  "mitochondria", "mito",
  "chloroplasts", "chloropl",
  "plastoglobuli", "plastoglob",
  "pyrenoid", "pyrenoid",
  "starch", "starch",
  "large_light_core_vesicle", "core_vesic",
  "food_vacuole", "food_vac",
  "lysosome", "lysosome",
  "reticulated_net", "retic_net",
  "crystal_rich_ret_net", "cryst_rich_RN",
  "electron_dense_sheets", "electr_sheets",
  "pusule", "pusule",
  "tubular_net", "tubul_net",
  "thrichocysts", "tricho",
  "eyespot", "eyespot",
  "rhabdosome", "rhabdo",
  "SER_whorls", "SER_whorls",
  "lipid_droplets", "lipid_drop",
  "fibrous_body", "fibrous_body",
  "peduncle", "peduncle",
  "flagellar_apparatus", "flagell_app",
  "silica_deposition_vesicle", "SDV",
  "putative_cell_division", "put_division",
  "symbiosis", "symbiosis",
  "undescribed_organelles", "undescribed",
  "beauty",  "beauty",
  "comments", "comments"
)

## columns describing the image; always taken from the TEM collection table
collection_identity_columns <- c(
  "name", "uri", "view", "grid", "site", "cell_id", "size_frac",
  "sampling_time", "source_name", "exclusive", "min_gray", "max_gray",
  "contrast_limits"
)

## columns with a fixed meaning in the personal tabs / main table
round_meta_columns <- c("round", "annotated_by", "validated_by")

## columns of old exports that are not annotation values
ignored_columns <- c("file", "run", "done_by")

## known misspellings of annotator names (lower case) -> correct name
name_aliases <- c(
  viktroria="viktoria",
  vitkoria="viktoria"
)

is_blank <- function(x) {
  is.na(x) | str_trim(as.character(x)) == ""
}

normalize_person <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  x[!is.na(x) & x == ""] <- NA_character_
  alias <- unname(name_aliases[x])
  if_else(!is.na(alias), alias, x)
}

## number of filled annotation fields per row
count_filled <- function(df) {
  cols <- intersect(annotation_columns$short, names(df))
  if (length(cols) == 0 || nrow(df) == 0) {
    return(integer(nrow(df)))
  }
  reduce(map(cols, ~as.integer(!is_blank(df[[.x]]))), `+`)
}

## annotation value columns of a table: everything that is not identity, meta or internal
value_columns <- function(df) {
  extra <- setdiff(
    names(df),
    c(collection_identity_columns, round_meta_columns, ignored_columns)
  ) %>%
    discard(~str_starts(.x, fixed(".")))
  c(annotation_columns$short, setdiff(extra, annotation_columns$short))
}

## standard column order of the main table
order_main_columns <- function(df) {
  df %>%
    select(
      any_of(collection_identity_columns), round, annotated_by,
      any_of(annotation_columns$short), everything()
    )
}

## standard column order of a personal tab
order_tab_columns <- function(df) {
  df %>%
    select(
      any_of(collection_identity_columns), round,
      any_of(setdiff(annotation_columns$short, "comments")),
      annotated_by, validated_by, any_of("comments"), everything()
    )
}

## Classify rows of personal tabs (live or from backups).
## Needs the columns name, annotated_by, validated_by (round optional) and
## .sheet / .run / .source (tab owner, timestamp, origin) added by the caller.
## .status: empty (nothing filled), in_progress (values but nobody signed), signed
classify_tab_rows <- function(df) {
  if (nrow(df) == 0) {
    return(tibble(
      name=character(), annotated_by=character(), validated_by=character(), round=character(),
      .sheet=character(), .run=character(), .source=character(), .n_filled=integer(),
      .ann=character(), .val=character(), .tab_round=integer(), .person=character(),
      .wanted_round=integer(), .status=character()
    ))
  }
  for (col in round_meta_columns) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_character_
    }
  }
  df$.n_filled <- count_filled(df)
  df %>%
    filter(!is_blank(name)) %>%
    mutate(
      .ann=normalize_person(annotated_by),
      .val=normalize_person(validated_by),
      .tab_round=suppressWarnings(as.integer(round)),
      ## the name in annotated_by / validated_by counts, not the tab it was entered in
      .person=case_when(
        .tab_round %in% 2L ~ .val,
        .tab_round %in% 1L ~ .ann,
        TRUE ~ coalesce(.val, .ann)
      ),
      .wanted_round=case_when(
        .tab_round %in% c(1L, 2L) ~ .tab_round,
        !is.na(.val) ~ 2L,
        TRUE ~ NA_integer_
      ),
      .status=case_when(
        .n_filled == 0 ~ "empty",
        is.na(.person) ~ "in_progress",
        TRUE ~ "signed"
      )
    )
}

## stable text signature of the annotation values (to compare versions)
value_signature <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if (length(cols) == 0) {
    return(rep("", nrow(df)))
  }
  df %>%
    select(all_of(cols)) %>%
    mutate(across(everything(), ~if_else(is_blank(.x), "", str_trim(as.character(.x))))) %>%
    pmap_chr(paste, sep="\u001f")
}

## Merge signed tab rows into the existing main table.
## existing:   main table (name, round, annotated_by, values), may be empty
## candidates: signed rows from classify_tab_rows()
## Rules:
##  - per image and person only one version is used: most filled fields, then the newest
##  - a person who already has a round of this image in main is not added again
##  - validated_by (or a tab round 2) claims round 2; otherwise the earliest person
##    takes round 1 and the next different person round 2
##  - a third person (both rounds taken) is reported, not added
## Returns list(rounds=<new main table>, issues=<tibble>)
resolve_rounds <- function(existing, candidates) {
  issues <- list()
  vcols <- unique(c(value_columns(existing), value_columns(candidates)))

  if (nrow(existing) == 0) {
    existing <- tibble(name=character(), round=integer(), annotated_by=character())
  }
  existing <- existing %>%
    mutate(
      round=as.integer(round),
      annotated_by=normalize_person(annotated_by)
    )

  if (nrow(candidates) == 0) {
    return(list(rounds=existing, issues=tibble()))
  }

  candidates <- candidates %>%
    mutate(.sig=value_signature(candidates, vcols))

  ## several differing versions of the same person
  multi <- candidates %>%
    group_by(name, .person) %>%
    filter(n_distinct(.sig) > 1) %>%
    summarise(
      detail=paste0(
        n(), " versions (", n_distinct(.sig), " differing), filled fields: ",
        paste(.n_filled, collapse="/"), "; sources: ", paste(unique(.source), collapse=", ")
      ),
      .groups="drop"
    )
  if (nrow(multi) > 0) {
    issues$multi <- multi %>%
      transmute(type="several_versions_kept_most_complete", name, person=.person, detail)
  }

  best <- candidates %>%
    group_by(name, .person) %>%
    mutate(
      .first_run=min(.run),
      .wanted_round=if (any(.wanted_round %in% 2L)) 2L else suppressWarnings(min(.wanted_round, na.rm=TRUE))
    ) %>%
    arrange(desc(.n_filled), desc(.run), .by_group=TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(.wanted_round=as.integer(replace(.wanted_round, is.infinite(.wanted_round), NA)))

  ## person already has a round of this image in main
  already <- best %>%
    inner_join(
      existing %>% mutate(.existing_sig=value_signature(existing, vcols)) %>%
        select(name, .person=annotated_by, .existing_round=round, .existing_sig),
      by=c("name", ".person")
    )
  if (nrow(already) > 0) {
    differing <- already %>% filter(.sig != .existing_sig)
    if (nrow(differing) > 0) {
      issues$already <- differing %>%
        transmute(
          type="already_in_main_tab_version_ignored", name, person=.person,
          detail=paste0("round ", .existing_round, " already in main; differing version from ", .source)
        )
    }
  }
  best <- best %>% anti_join(already %>% select(name, .person), by=c("name", ".person"))

  added <- list()
  conflicts <- list()
  for (img in unique(best$name)) {
    taken <- existing %>% filter(name == img)
    r1 <- taken$annotated_by[taken$round == 1L]
    r2 <- taken$annotated_by[taken$round == 2L]
    cand <- best %>%
      filter(name == img) %>%
      arrange(desc(.wanted_round %in% 2L), .first_run)
    for (i in seq_len(nrow(cand))) {
      row <- cand[i, , drop=FALSE]
      person <- row$.person
      target <- NA_integer_
      if (row$.wanted_round %in% 2L) {
        if (length(r2) == 0 && !person %in% r1) target <- 2L
      } else if (length(r1) == 0 && !person %in% r2) {
        target <- 1L
      } else if (length(r2) == 0 && !person %in% r1) {
        target <- 2L
      }
      if (is.na(target)) {
        conflicts[[length(conflicts) + 1]] <- tibble(
          type="both_rounds_taken_not_added", name=img, person=person,
          detail=paste0(
            "round 1: ", paste(r1, collapse=","), ", round 2: ", paste(r2, collapse=","),
            "; source: ", row$.source
          )
        )
        next
      }
      if (target == 1L) r1 <- person else r2 <- person
      added[[length(added) + 1]] <- row %>% mutate(round=target, annotated_by=person)
    }
  }
  if (length(conflicts) > 0) {
    issues$conflicts <- bind_rows(conflicts)
  }

  new_rows <- bind_rows(added)
  if (nrow(new_rows) > 0) {
    new_rows <- new_rows %>%
      select(any_of(collection_identity_columns), round, annotated_by, any_of(vcols), .source, .run)
  }

  missing_r1 <- bind_rows(existing, new_rows) %>%
    group_by(name) %>%
    filter(any(round == 2L), !any(round == 1L)) %>%
    ungroup()
  if (nrow(missing_r1) > 0) {
    issues$missing_r1 <- missing_r1 %>%
      transmute(type="round2_without_round1", name, person=annotated_by, detail="validation found but no annotation round")
  }

  rounds <- bind_rows(existing, new_rows %>% select(-any_of(c(".source", ".run"))))
  list(rounds=rounds, issues=bind_rows(issues), added=new_rows)
}

## replace identity columns of the main table with the values from the collection table;
## rows of images that are not in the collection table keep their old values
attach_identity <- function(rounds, collection) {
  if (nrow(rounds) == 0) {
    return(rounds)
  }
  id_cols <- intersect(setdiff(collection_identity_columns, "name"), names(collection))
  in_coll <- rounds %>%
    filter(name %in% collection$name) %>%
    select(-any_of(id_cols)) %>%
    left_join(collection %>% select(name, all_of(id_cols)), by="name")
  not_in_coll <- rounds %>% filter(!name %in% collection$name)
  bind_rows(in_coll, not_in_coll)
}

## one row per image of the collection table: who did round 1 and round 2
make_status_table <- function(collection, rounds) {
  people_by_round <- function(r) {
    rounds %>%
      filter(round == r) %>%
      distinct(name, .keep_all=TRUE) %>%
      select(name, person=annotated_by)
  }
  collection %>%
    select(any_of(collection_identity_columns)) %>%
    left_join(people_by_round(1L) %>% rename(annotated_by=person), by="name") %>%
    left_join(people_by_round(2L) %>% rename(validated_by=person), by="name") %>%
    arrange(site, name)
}

sort_rounds <- function(rounds) {
  if (!"site" %in% names(rounds)) {
    rounds$site <- rep(NA_character_, nrow(rounds))
  }
  rounds %>%
    mutate(round=as.integer(round)) %>%
    arrange(site, name, round) %>%
    order_main_columns()
}
