


library(tidyverse)

trec24 <- read_csv("C:/projects/trec_biobank/metab/TREC_Plankton_wbrMetaB_wbrTRECAML_wbr40um_wbr2024_dataset_files_export.csv")

trec23 <- read_csv("C:/projects/trec_biobank/metab/TREC_Plankton_wbrMetaB_wbrTRECAML_wbr40um_dataset_files_export.csv")

vigo <- read_csv("C:/projects/trec_biobank/metab/TREC_CELLBLOOM_wbrVigo_wbrmetabarcoding_dataset_files_export.csv")


neg_control <- c("SAMEA122159368",
  "SAMEA122248694", "SAMEA122248734", "SAMEA122248735", "SAMEA122248737" ,"SAMEA122248738")


full_dataset_table <- bind_rows(trec24, trec23, vigo) %>%
  mutate(type=ifelse(Samples %in% neg_control, "control", "environmental")) %>%
  select(-Type, -Format)

write_tsv(full_dataset_table, file="C:/projects/trec_biobank/metab/sampleinfo.tsv")

test %>% filter(Samples %in% neg_control) %>% View()





library(ijtiff)
img <- read_tif("A:/Marco/Chandni_mitodino/CellE_chromosomes.tif", drop = TRUE)

img <- readImage("A:/Marco/Chandni_mitodino/CellE_chromosomes.tif", all=T)



googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1Wq9Fxh_tiPQn20F8sxrXquyDgU_ZkhP248v6IvHcdrY/edit?gid=0#gid=0")


img[1,1,1]


array_img <- lapply(seq(1, length(img)), function(n){
  print(n)
  return(img[[n]][,,2])
  
})


test <- array_img[1:3]


t <- array(array_img, dim=c(1894, 1836, 1188))






df_asv_per_sample <- read_tsv("A:/marco/projects/trec_metaB/prepared_input/asvs_per_sample.tsv")
df_asv_taxonomy <- read_tsv("A:/marco/projects/trec_metaB/prepared_input/asv_taxonomy.tsv")
df_sampling_sites <- read_tsv("A:/marco/projects/trec_metaB/prepared_input/sampling_sites.tsv")






taxonomy <- df_asv_taxonomy %>%
  filter(str_detect(level_4, "inoflag")) %>%
  select(1,3:11) 


level_names <- c(
  level_1 = "Domain",
  level_2 = "Supergroup",
  level_3 = "Division",
  level_4 = "Phylum",
  level_5 = "Class",
  level_6 = "Order",
  level_7 = "Family",
  level_8 = "Genus",
  level_9 = "Species"
)

missing_taxonomy <- taxonomy %>%
  summarise(
    across(
      starts_with("level_"),
      ~ sum(is.na(.) | . == "NA")
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "level",
    values_to = "n_missing"
  ) %>%
  mutate(
    taxonomic_rank = recode(level, !!!level_names),
    total_rows = nrow(taxonomy),
    percent_missing = 100 * n_missing / total_rows
  ) %>%
  select(
    level,
    taxonomic_rank,
    n_missing,
    total_rows,
    percent_missing
  )%>%
  mutate(
    taxonomic_rank = factor(
      taxonomic_rank,
      levels = c(
        "Domain", "Supergroup", "Division", "Phylum",
        "Class", "Order", "Family", "Genus", "Species"
      )
    )
  )

trec_p <- ggplot(missing_taxonomy %>% filter(percent_missing!=0),
       aes(x=taxonomic_rank, y=percent_missing))+
  geom_col()+
  theme_bw()

png(
  filename = "C:/projects/trec_biobank/analysis/metab/rstudio/plots/missing_taxonomy.png",
  width = 3,
  height = 3,
  units = "in",
  res = 300,
  bg = "white"
)

trec_p

dev.off()





%>%
  mutate(
    level = factor(
      level,
      levels = c(
        "Domain", "Supergroup", "Division", "Phylum",
        "Class", "Order", "Family", "Genus", "Species"
      )
    )
  ) %>%
  ggplot(aes(x = level, y = percent_missing)) +
  geom_col(fill = "#4E79A7", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent_missing)),
    vjust = -0.3,
    size = 4
  ) +
  labs(
    x = "Taxonomic rank",
    y = "Missing annotations (%)"
  ) +
  ylim(0, max(missing_taxonomy$percent_missing) * 1.1) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )














fac <- c("HR", "AS", "AC", "PL", "P", "GI", "HeLa", "scrambled")




input %>%
  select(-1) %>%
  pivot_longer(cols = c("HR", "AS", "AC", "PL", "P", "GI", "HeLa"), names_to="species") %>%
  mutate(species=factor(species, levels=fac),
         probes=factor(`probes/species`, levels=fac),
         ) %>%
  ggplot(.,
         aes(x=species, y=probes, fill=as.character(value)))+
  geom_tile()+
  scale_fill_manual(values=c("orange", "lightblue", "blue"))










input <- googlesheets4::read_sheet(
  "https://docs.google.com/spreadsheets/d/1Wq9Fxh_tiPQn20F8sxrXquyDgU_ZkhP248v6IvHcdrY/edit?gid=0#gid=0",
  sheet="optimal")




fac <- c("HR", "AS", "AC", "PL", "P", "GI", "HeLa", "scrambled")

species_labels <- c(
  HR   = "Heterocapsa rotundata",
  AS   = "Akashiwo sanguinea",
  AC   = "Amphidinium carterae",
  PL   = "Pyrocystis lunula",
  P    = "Protodinium sp.",
  GI   = "Gymnodinium impudicum",
  HeLa = "HeLa",
  scrambled = "Scrambled"
)

plot <- input %>%
  select(-1) %>%
  pivot_longer(
    cols = c("HR", "AS", "AC", "PL", "P", "GI", "HeLa"),
    names_to = "species"
  ) %>%
  mutate(
    species = factor(species, levels = fac),
    probes  = factor(`probes/species`, levels = fac)
  ) %>%
  ggplot(aes(x = species, y = probes, fill = as.character(value))) +
  geom_tile() +
  scale_fill_manual(values = c("orange", "lightblue", "blue")) +
  scale_x_discrete(labels = species_labels) +
  scale_y_discrete(labels = species_labels) + 
  scale_fill_manual(
    name = "Signal",
    values = c(
      "0" = "#E0E0E0",
      "1" = "#FDB863",
      "2" = "#1B7837"
    ),
    labels = c(
      "0" = "No signal",
      "1" = "Weak signal",
      "2" = "Strong signal"
    )
  )+
  theme(
    axis.text.x = element_text(angle=315, hjust=0)
  )



png(
  filename = "C:/projects/trec_biobank/analysis/metab/rstudio/plots/probe_signal_heatmap_optimal.png",
  width = 1950,
  height = 1600,
  res = 300
)
plot
dev.off()





























