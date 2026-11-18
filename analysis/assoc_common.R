## Shared by the association analyses (04/05 organelles x organelles, 07/08 organelles x metadata)
## and 06_organelle_clustering.R: feature selection, tests and plot helpers

## ---- feature selection (used by 04, 05 and 06) ----
## organelles to test / cluster on: remove / comment out what should not be used
## (features that are hard to judge by the study design stay out -> fewer tests, fewer false positives)
use_features <- c(
  #"cell_cover", ## noisy column... a lot of different things count as cell cover 
  #"nucleus",    ## we know that every cell has it
  #"nucleolus",  ## unlikely to cut
  # "pusule", ## unlikely to cut
  #"flagell_app", ## unlikely to cut
  #"symbiosis",  ## extremely rare so probably not suited for this analysis but actually interesting to correlate with organelles
  #"eyespot", ## unlikely to cut
  #"put_division", 
  # "SDV", 
  #"ER", "mito", "golgi", ## every cell has it
  ## to be discussed: 
  # "lysosome", ## a bit difficult to annotate... noidsy ... per person bias
  # "fibrous_body", some cells have a lot i remeber so you would always cut it but some also not?
  "pyrenoid",
  "tricho", 
  "chloropl", 
  "plastoglob", 
  "retic_net", 
  "electr_sheets", 
  "tubul_net",
  "food_vac", 
  "starch", 
  "cryst_rich_RN", 
  "rhabdo", 
  "core_vesic",  ## not sure?
  "lipid_drop", ## difficult to annotae, so noisy but could tell us function about energy storage
  "SER_whorls" ## not sure?
  )
include_undescribed <- TRUE    # also use the undescribed_<n> structures
min_feature_cells   <- 10      # features present in fewer cells are dropped before all analyses (04-09)

## nominal columns (some coded as numbers, but without order) -> treat as categorical
nominal_cols <- c("site", "size_frac", "sampling_time", "life_status", "major_group")

## organelle features of the analyses: use_features (+ undescribed_<n> if include_undescribed),
## without the features present in fewer than min_feature_cells cells (too rare to trust)
assoc_features <- function(dino, verbose = FALSE) {
  unknown <- setdiff(use_features, names(dino))
  if (length(unknown) > 0) stop("use_features not in the data: ", paste(unknown, collapse = ", "))
  feats <- c(use_features, if (include_undescribed) undescribed_cols(dino))
  n_present <- map_int(feats, ~ sum(suppressWarnings(as.numeric(dino[[.x]])) == 1, na.rm = TRUE))
  rare <- n_present < min_feature_cells
  if (verbose && any(rare))
    message(sprintf("%d features in < %d cells dropped: %s", sum(rare), min_feature_cells,
                    paste0(feats[rare], " (", n_present[rare], ")", collapse = ", ")))
  feats[!rare]
}

## metadata the organelles are tested against (07/08); columns with only one value are dropped
assoc_metadata <- function(cor_data) {
  keep(intersect(meta_feats, names(cor_data)), ~ n_distinct(na.omit(cor_data[[.x]])) > 1)
}

## metadata (meta_feats) + organelle features
build_cor_data <- function(dino) {
  dino %>%
    select(any_of(meta_feats), all_of(assoc_features(dino))) %>%
    mutate(major_group = na_if(as.character(major_group), "?"),
           across(any_of(nominal_cols), as.character))
}

is_binary <- function(x) is.numeric(x) && all(na.omit(x) %in% c(0, 1))

## ---- tests (04 organelles x organelles, 07 organelles x metadata) ----
## test by column type, NAs dropped per pair:
##   0/1 vs 0/1           fisher 2x2, effect = phi, + conditional MLE odds ratio with 95% CI
##   numeric vs numeric   spearman (beauty is ordinal)
##   any categorical      fisher RxC with simulated p, effect = Cramer's V
test_pair <- function(var1, var2, data) {
  d <- data %>% select(x = all_of(var1), y = all_of(var2)) %>% drop_na()
  out <- tibble(var1, var2, n = nrow(d),
                n1 = nrow(d %>% filter(x == 1)),
                n2 = nrow(d %>% filter(y == 1)),
                test = NA_character_, estimate = NA_real_,
                var1_present = NA_character_, var1_absent = NA_character_,
                odds_ratio = NA_real_, or_ci_low = NA_real_, or_ci_high = NA_real_, p_value = NA_real_)
  if (n_distinct(d$x) < 2 || n_distinct(d$y) < 2) return(out)

  if (is_binary(d$x) && is_binary(d$y)) {
    ## 0/1 vs 0/1: exact test on the 2x2 table, effect = phi coefficient (-1..1)
    ## + conditional MLE odds ratio with 95% CI (only defined for 2x2)
    ## counts behind the OR as "var2 present / var2 absent", split by var1 present / absent
    ## OR ~ (a/b) / (c/d)
    tab <- table(factor(d$x, levels = 0:1), factor(d$y, levels = 0:1))
    res <- fisher.test(tab)
    out %>% mutate(test = "fisher 2x2", estimate = cor(d$x, d$y),
                   var1_present = paste0(tab["1", "1"], "/", tab["0", "0"]),
                   var1_absent  = paste0(tab["0", "1"], "/", tab["1", "0"]),
                   x1y1 = tab["1", "1"], x0y0 = tab["0", "0"], x0y1 = tab["0", "1"], x1y0 = tab["1", "0"],
                   odds_ratio = unname(res$estimate), or_ci_low = res$conf.int[1], or_ci_high = res$conf.int[2],
                   p_value = res$p.value)
  } else if (is.numeric(d$x) && is.numeric(d$y)) {
    ## ordinal (beauty) vs numeric: rank correlation (-1..1)
    res <- cor.test(d$x, d$y, method = "spearman", exact = FALSE)
    out %>% mutate(test = "spearman", estimate = unname(res$estimate), p_value = res$p.value)
  } else {
    ## at least one categorical: fisher with simulated p, effect = Cramer's V (0..1)
    tab <- table(d$x, d$y)
    res <- fisher.test(tab, simulate.p.value = TRUE, B = 1e4)
    chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
    out %>% mutate(test = "fisher RxC", estimate = unname(sqrt(chi2 / (sum(tab) * (min(dim(tab)) - 1)))),
                   p_value = res$p.value)
  }
}

## all pairs (tibble var1, var2), Bonferroni over the pairs of this analysis only
run_assoc_tests <- function(pairs, data) {
  set.seed(42)
  pmap_dfr(pairs, function(var1, var2) test_pair(var1, var2, data)) %>%
    mutate(p_adj = p.adjust(p_value, method = "bonferroni")) %>%
    arrange(p_adj)
}

## writes <prefix>_results.tsv/.rds and reports the numbers
save_assoc_results <- function(results, prefix, table_dir) {
  write_tsv(results, file.path(table_dir, paste0(prefix, "_results.tsv")))
  saveRDS(results, file.path(table_dir, paste0(prefix, "_results.rds")))
  message(sprintf("%s: %d pairs tested, %d significant (Bonferroni p < 0.05)", prefix,
                  sum(!is.na(results$p_value)), sum(results$p_adj < 0.05, na.rm = TRUE)))
}

read_assoc_results <- function(prefix, table_dir, step) {
  file <- file.path(table_dir, paste0(prefix, "_results.rds"))
  if (!file.exists(file)) stop("Run ", step, " first (missing ", file, ")")
  readRDS(file)
}

## ---- plot helpers (05, 08) ----
stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
est_names <- c("fisher 2x2" = "phi", "spearman" = "rho", "fisher RxC" = "Cramer's V")

feature_group <- function(f) {
  factor(case_when(f %in% chloro_feats ~ "chloroplast",
                   f %in% common_feats ~ "common",
                   f %in% meta_feats   ~ "metadata",
                   str_detect(f, "^undescribed_") ~ "undescribed",
                   TRUE ~ "specific"),
         c("chloroplast", "common", "specific", "metadata", "undescribed"))
}

## composition of var2 within the categories of var1 (variable with more levels on the x axis)
plot_cat_pair <- function(var1, var2, p_adj, estimate, test, data) {
  d <- data %>% select(x = all_of(var1), y = all_of(var2)) %>% drop_na() %>%
    mutate(across(c(x, y), as.factor))
  if (nlevels(d$y) > nlevels(d$x)) {
    d <- rename(d, x = y, y = x)
    tmp <- var1; var1 <- var2; var2 <- tmp
  }
  ggplot(d, aes(x, fill = y)) +
    geom_bar(position = "fill") +
    geom_text(data = count(d, x), aes(x, y = 1.02, label = n), inherit.aes = FALSE, size = 3, vjust = 0) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = "Set2", name = var2) +
    labs(x = var1, y = "fraction of cells", title = paste(var1, "vs", var2),
         subtitle = sprintf("%s = %.2f, p_adj = %.2g (n per bar on top)", est_names[test], estimate, p_adj)) +
    theme_minimal()
}

## all composition plots in one multi-page pdf, one png per pair
save_cat_plots <- function(plots, prefix, plot_dir) {
  if (length(plots) == 0) return(invisible(NULL))
  pdf(file.path(plot_dir, paste0(prefix, ".pdf")), width = 7, height = 5)
  walk(plots, print)
  dev.off()
  iwalk(plots, ~ ggsave(file.path(plot_dir, sprintf("%s_%02d.png", prefix, .y)), .x,
                        width = 7, height = 5, dpi = 300, bg = "white"))
}

## contingency table of every significant pair, ordered by p_adj, as one faceted plot
save_contingency_tables <- function(sig, data, name, plot_dir) {
  if (nrow(sig) == 0) return(invisible(NULL))
  ct_long <- sig %>%
    select(var1, var2, test, estimate, p_adj) %>%
    pmap_dfr(function(var1, var2, test, estimate, p_adj) {
      data %>%
        select(x = all_of(var1), y = all_of(var2)) %>%
        drop_na() %>%
        mutate(across(c(x, y), as.character)) %>%
        count(x, y) %>%
        complete(x, y, fill = list(n = 0L)) %>%
        mutate(frac_x = n / sum(n), .by = x) %>%
        mutate(pair = sprintf("x: %s | y: %s\n%s = %.2f, p_adj = %.1e", var1, var2, est_names[test], estimate, p_adj))
    }) %>%
    mutate(pair = factor(pair, unique(pair)))

  n_col <- ceiling(sqrt(n_distinct(ct_long$pair)))
  n_row <- ceiling(n_distinct(ct_long$pair) / n_col)

  p_tables <- ggplot(ct_long, aes(x, y, fill = frac_x)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = n), size = 3) +
    facet_wrap(~ pair, scales = "free", ncol = n_col) +
    scale_fill_gradient(low = "white", high = "#d95f02", limits = c(0, 1), labels = scales::percent,
                        name = "share within\nx category") +
    labs(x = NULL, y = NULL, title = "Contingency tables of all significant pairs",
         subtitle = "numbers = cell counts, colour = distribution of y within each x category (column)") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), strip.text = element_text(size = 7, hjust = 0))

  save_plot(p_tables, name, 3 * n_col, 2.6 * n_row + 1, plot_dir)
}

## ---- MoBIE collection tables (06 clusters, 09 negative pairs; overview figures in step 10) ----
## columns of a collection table (MoBIE + image metadata); view/grid are set to the group
collection_cols <- c("name", "uri", "view", "grid", "site", "cell_id", "size_frac", "sampling_time",
                     "source_name", "exclusive", "min_gray", "max_gray", "contrast_limits")

## spreadsheet for the collection tables (--collection_sheet_url / env TEM_COLLECTION_SHEET_URL);
## empty = only the local tsv files
collection_sheet <- function(opt) {
  opt$collection_sheet_url %||%
    Sys.getenv("TEM_COLLECTION_SHEET_URL",
               "https://docs.google.com/spreadsheets/d/18b2AVslU7sOQOTRRIwxGHKPtrh00ezT-byOrQpiVo8Y/edit?gid=0#gid=0")
}

## write named tables as tabs (name = tab) into the collection spreadsheet
push_collection_tabs <- function(tables, opt) {
  url <- collection_sheet(opt)
  if (is.na(url) || url == "") return(invisible(NULL))
  suppressPackageStartupMessages(library(googlesheets4))
  if (!file.exists(opt$google_key)) stop("Google key not found: ", opt$google_key)
  gs4_auth(path = opt$google_key)
  tryCatch({
    iwalk(tables, ~ write_sheet(.x, ss = url, sheet = .y))
    message("Written to ", url, ": ", paste(names(tables), collapse = ", "))
  }, error = function(e) {
    stop("Writing to the collection spreadsheet failed (the service account ",
         jsonlite::read_json(opt$google_key)$client_email,
         " needs edit access to it): ", conditionMessage(e))
  })
}
