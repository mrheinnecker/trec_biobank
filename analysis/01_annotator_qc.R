## Annotator QC
## 1) beauty score distribution per annotator
## 2) uncertain entries ("?") that were counted as absent: table + overview figure
## 3) agreement between the first annotation (round 1) and the validation (round 2);
##    only when round 2 is selected (--rounds 1,2)
## annotated_by = the person who did that round

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- NA_character_   # NA = read the main sheet from Google Sheets; or a tsv snapshot
default_outdir     <- "C:/projects/tem_classification"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/tem_classification/analysis"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)

all_annotations <- load_annotations(opt)

## ---- beauty per annotator ----
beauty_plot <- all_annotations %>%
  filter(!is.na(beauty)) %>%
  count(annotated_by, beauty) %>%
  ggplot(aes(x = annotated_by, y = n, fill = beauty)) +
  geom_col(position = "fill") +
  scale_fill_gradientn(colors = c("black", "blue", "deeppink"))

save_plot(beauty_plot, "beauty", 5, 4, opt$plot_dir)

## ---- uncertain entries ----
## every "?" entry that was set to absent (see clean_uncertain() in common.R)
uncertain <- attr(all_annotations, "uncertain")
write_tsv(uncertain, file.path(opt$table_dir, "uncertain_entries.tsv"), na = "")

if (nrow(uncertain) > 0) {
  uncertain_counts <- uncertain %>%
    count(column, annotated_by) %>%
    mutate(column = fct_reorder(column, n, .fun = sum))

  p_uncertain <- ggplot(uncertain_counts, aes(x = annotated_by, y = column, fill = n)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = n, colour = n > 0.6 * max(n)), size = 3.5) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1a1a1a"), guide = "none") +
    scale_fill_gradient(low = "#cde2fb", high = "#184f95", name = "entries") +
    labs(x = "annotated by", y = NULL,
         title = "Uncertain entries counted as absent",
         subtitle = sprintf("%d entries in %d cells", nrow(uncertain), n_distinct(uncertain$name)),
         caption = "details: tables/uncertain_entries.tsv") +
    theme_minimal() +
    theme(panel.grid = element_blank(), plot.title.position = "plot")

  save_plot(p_uncertain, "uncertain_entries", 3.5 + 0.9 * n_distinct(uncertain_counts$annotated_by),
            1.5 + 0.35 * n_distinct(uncertain_counts$column), opt$plot_dir)
}

## ---- annotation vs validation ----
if (!2 %in% opt$rounds) {
  message("Round 2 not selected (--rounds): annotation vs validation agreement skipped")
} else {
  feature_cols <- all_annotations %>%
    select(life_status:last_col()) %>%
    select(-any_of(c("comments", "annotated_by", "validated_by", "round", "n_cells"))) %>%
    names()

  ## cells with both rounds: one row per cell and feature, value of round 1 vs round 2
  difference <- all_annotations %>%
    group_by(name) %>%
    filter(all(c(1, 2) %in% round)) %>%
    ungroup() %>%
    select(name, round, annotated_by, all_of(feature_cols)) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(all_of(feature_cols), names_to = "key", values_to = "value") %>%
    pivot_wider(names_from = round, values_from = c(value, annotated_by)) %>%
    rename(round_1 = value_1, round_2 = value_2, annotator = annotated_by_1, validator = annotated_by_2) %>%
    mutate(same = round_1 == round_2)

  ## features ordered by number of disagreements
  lev <- difference %>%
    mutate(dff = !same) %>%
    group_by(key) %>%
    summarise(n = sum(dff, na.rm = TRUE)) %>%
    arrange(desc(n)) %>%
    pull(key)

  p_acceptance <- ggplot(difference %>% mutate(key = factor(key, levels = lev)),
                         aes(y = name, x = key, fill = same)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 315, hjust = 0))

  save_plot(p_acceptance, "acceptance", 13, 4, opt$plot_dir)
  write_tsv(difference, file.path(opt$table_dir, "validation_agreement.tsv"))
}
