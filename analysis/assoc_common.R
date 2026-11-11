## Shared by 04_associations.R (tests) and 05_association_plots.R (plots):
## the variable set of the association analysis

## nominal columns (some coded as numbers, but without order) -> treat as categorical
nominal_cols <- c("site", "size_frac", "sampling_time", "life_status", "major_group")

## metadata + all annotated features (incl. undescribed_<n>), without bookkeeping / free text columns
build_cor_data <- function(dino) {
  dino %>%
    select(site, size_frac, sampling_time, life_status:last_col()) %>%
    select(-comments, -annotated_by, -validated_by, -taxo_class, -undescribed, -file, -run, -done_by) %>%
    mutate(major_group = na_if(as.character(major_group), "?"),
           across(any_of(nominal_cols), as.character))
}

is_binary <- function(x) is.numeric(x) && all(na.omit(x) %in% c(0, 1))
