library(tidyverse)
library(googlesheets4)
library(googledrive)

google_key <- "/g/schwab/marco/repos/trec_biobank/trec-tem-screen-e98a2e03f58b.json"
gs4_auth(path=google_key)
drive_auth(path=google_key)
tem_collection_table_url <- "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1643802951#gid=1643802951"

## load collection tbale: this one has all datasets that were correctly converted and are sitting on the s3
tem_collection_table <- read_sheet(tem_collection_table_url, sheet="tem_collection_table") %>%
  select(1:13)

## load old annotations from viktoria for now (until final adaptaion is done)
image_log_url <- "https://docs.google.com/spreadsheets/d/143uVeeJ72SQE5eK01lzWYCEiT7pJUF3lX7hJl3R9s9I/edit?gid=2107269356#gid=2107269356"

image_log <- read_sheet(image_log_url, sheet="image_log") 

## filter out some that we need to correct before annotations:
## These are either those that are at the moment done by Viktoria (here we need to adapt them to the final column naming conventions)

annotated_by_viktoria <- image_log %>%
  ## an annoatted column works here
  filter(!is.na(nucleus))

multiple_cells <- image_log %>%
  filter(cell_count>1)

remove_for_now <- c(annotated_by_viktoria$shortname, multiple_cells$shortname) %>% unique()


## now i create a dataframe to stash the removed ones until manula correczion is done

annotation_main <- tem_collection_table %>%
  filter(!name %in% remove_for_now)

to_be_reannoatted_viktoria <- tem_collection_table %>%
  filter(name %in% annotated_by_viktoria$shortname) %>%
  left_join(image_log %>%
              select(name=shortname,
                     
                     n_cells=cell_count, 
                     life_status, 
                     cell_cover=theca, 
                     nucleus, 
                     nucleolus, 
                     golgi=golgi_apparatus, 
                     mito=mitochondria, 
                     chloropl=chloroplasts, 
                     plastoglob=plastoglobuli, 
                     starch, 
                     core_vesic=vacuole, 
                     ER, 
                     retic_net=`reticulated_network/crystalline_compartment`, 
                     tubul_net=tubular_network, 
                     pusule, 
                     electr_sheets=electron_dense_sheets, 
                     rhabdo=rhabdosome, 
                     eyespot, 
                     flagell_app=`basal_body/flagellum`, 
                     symbiosis=`symbiotic/parasitism`, tmp_nin_ident=`non-identified`, 
                     comments=comment
                     
                     ), 
            by="name")



#Zu den kleinen Änderungen in der Tabelle, ich hab thecal_struc zwischen cell_cover und nucleus 
#Und peduncle zwischen fibrous_body und flagell_app
#Und ich hab lysosome hinter food_vac verschoben


annotation_columns <- tribble(
  ~full, ~short,
  "cell_count",  "n_cells",
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
  "annotated_by", "annotated_by",
  "validated_by", "validated_by",
  "comments", "comments"
  
)  


# viktorias suggestions for column order
# "
# n_cells
# life_status
# major_group
# taxo_class
# cell_cover
# nucleus
# nucleolus
# golgi
# ER
# mito
# chloropl
# plastoglob
# pyrenoid
# starch 
# core_vesic
# food_vac
# retic_net
# cryst_rich_RN
# electr_sheets
# pusule
# tubul_net
# tricho
# 
# eyespot
# rhabdo
# lysosome
# SER_whorls
# lipid_drop
# fibrous_body
# flagell_app
# SDV
# put_division
# symbiosis
# undescribed
# 
# 
# "
# 



all_cols_df <- tibble(
  !!!setNames(rep(list(NA), nrow(annotation_columns)),
              annotation_columns$short)
)


emtpy_full <- bind_cols(annotation_main, all_cols_df)
# 
# 
# main_annotations_url <- "https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=0#gid=0"
# 
# write_sheet(emtpy_full, ss=main_annotations_url, sheet="main")
# 
# 

filt_anno_cols <- annotation_columns %>%
  filter(!short %in% names(to_be_reannoatted_viktoria))

all_cols_df_only_those_toberedone <- tibble(
  !!!setNames(rep(list(NA), nrow(filt_anno_cols)),
              filt_anno_cols$short)
)


to_be_redone_by_viktoria <- bind_cols(to_be_reannoatted_viktoria, all_cols_df_only_those_toberedone)

ordered_redone <- to_be_redone_by_viktoria[, c("name", colnames(emtpy_full %>% select(-name)))]

# 
# main_annotations_url <- "https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=694322446#gid=694322446"
# 
# write_sheet(ordered_redone, ss=main_annotations_url, sheet="viktoria")
# 
# 


## now i also need to add them to the main table

to_be_added_main_tab <- tem_collection_table %>%
  filter(name %in% annotated_by_viktoria$shortname) %>%
  mutate(exclusive=as.character(exclusive))



main_annotations_url <- "https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=694322446#gid=694322446"


main <- read_sheet(main_annotations_url, sheet="main")

# 
# combined <- bind_rows(main, to_be_added_main_tab)
# 
# write_sheet(combined, ss=main_annotations_url, sheet="main")


