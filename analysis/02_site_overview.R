## Site overview: number of annotated cells per site, split by major group
## Uses the first annotation of every cell (all major groups).

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

annotations_first <- load_annotations(opt$input) %>% dedup_first_run()

tot_cells_per_site <- annotations_first %>%
  group_by(site) %>%
  tally()

site_major_group <- annotations_first %>%
  left_join(group_lookup, by = c("major_group" = "code")) %>%
  group_by(site, group) %>%
  tally()

p_ov_site_major_group <- ggplot(site_major_group, aes(x = site, y = n, fill = group)) +
  geom_col() +
  geom_text(data = tot_cells_per_site, inherit.aes = FALSE, aes(x = site, y = n, label = n), vjust = 0) +
  theme_bw()

save_plot(p_ov_site_major_group, "ov_site_major_group", 5, 3, opt$plot_dir)
write_tsv(site_major_group, file.path(opt$table_dir, "site_major_group_counts.tsv"))
