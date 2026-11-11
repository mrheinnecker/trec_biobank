## Site overview: number of annotated cells per site, split by major group
## Uses the first annotation of every cell (all major groups).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
opt <- parse_args()

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
