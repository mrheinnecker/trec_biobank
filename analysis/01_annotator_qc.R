## Annotator QC
## 1) beauty score distribution per annotator
## 2) agreement between the first annotation and the validation of validated cells
## Uses all annotation rows (no dedup), since it compares annotation rounds.

## ---- paths ----
## hard-coded for running line by line (RStudio); when run with Rscript,
## --input / --outdir (or the env vars TEM_ANNOT_INPUT / TEM_ANNOT_OUTDIR) take precedence
default_input      <- "C:/projects/tem_classification/all_annotations_2026-10-15.tsv"
default_outdir     <- "C:/projects/tem_classification"   # -> plots/ and tables/ in here
default_script_dir <- "C:/repos/tem_classification/analysis"

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir  <- if (is.na(script_file)) default_script_dir else dirname(normalizePath(sub("^--file=", "", script_file)))
source(file.path(script_dir, "common.R"))
opt <- parse_args(default_input = default_input, default_outdir = default_outdir)

all_annotations <- load_annotations(opt$input)

## validator if validated, annotator otherwise
add_done_by <- function(annot) {
  annot %>%
    mutate(done_by = case_when(
      !is.na(validated_by) ~ validated_by,
      !is.na(annotated_by) ~ annotated_by
    ))
}

## ---- beauty per annotator ----
beauty_plot <- all_annotations %>%
  filter(!is.na(beauty)) %>%
  add_done_by() %>%
  group_by(done_by, beauty) %>%
  tally() %>%
  mutate(
    done_by = factor(done_by, levels = c("chandni", "marco", "Karel", "viktoria", "yannick", "viktroria", "vitkoria"))
  ) %>%
  ggplot(aes(x = done_by, y = n, fill = beauty)) +
  geom_col(position = "fill") +
  scale_fill_gradientn(colors = c("black", "blue", "deeppink"))

save_plot(beauty_plot, "beauty", 5, 4, opt$plot_dir)

## ---- annotation vs validation ----
validated <- all_annotations %>%
  filter(!is.na(validated_by)) %>%
  pull(source_name)

difference <- all_annotations %>%
  filter(source_name %in% validated, !is.na(beauty)) %>%
  add_done_by() %>%
  select(name, done_by, 14:55) %>%
  select(-run, -comments, -file, -annotated_by, -validated_by) %>%
  mutate_all(.funs = as.character) %>%
  pivot_longer(cols = names(select(., 3:38)), names_to = "key", values_to = "value") %>%
  pivot_wider(names_from = "done_by", values_from = "value") %>%
  mutate(same = marco == viktoria)

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
