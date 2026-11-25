library(grid)
library(tidyverse)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(grid)  # for grobs
library(gridExtra)  # just in case
library(maps)

source("/g/schwab/marco/repos/trec_metaB/fncts.R")


input_tables <- load_required_tables()

p_list=lapply(sort(input_tables$site_mapping$site %>% unique()), 
              get_stacked_barplot_per_site, 
              df_asvs=input_tables$df_asvs, 
              site_mapping=input_tables$site_mapping,
              mapping=input_tables$mapping)



cities <- 
    tibble(
      name = c("Roscoff", "Kristineberg", "Villefranche",
               "Tallinn", 
               "Porto", 
               #"Barcelona", 
               #"Naples", "Athens", 
               "Reykjavik", "Bilbao"),
      lon = c(-3.9833, 11.4394, 7.3120,
              24.7536, -8.6110, #2.1734,
              #14.2681, 23.7275, 
              -21.8277, -2.9349),
      lat = c(48.7269, 58.2375, 43.7040,
              59.4370, 41.1496, #41.3851,
              #40.8522, 37.9838, 
              64.1283, 43.2630)
    )   %>%
  mutate(name = factor(name, levels = levels(input_tables$site_mapping$site))) %>%
  arrange(name)


ortho_crs <- "+proj=ortho +lat_0=50 +lon_0=10"
cities_sf <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)
cities_proj <- st_transform(cities_sf, crs = ortho_crs)

# Extract projected coordinates
coords <- st_coordinates(cities_proj)
cities_coords <- cbind(cities, coords)  # Add X/Y to original data

# Width and height in projection units (adjust as needed)
world <- ne_countries(scale = "medium", returnclass = "sf")
world_proj <- st_transform(world, crs = ortho_crs)
# Orthographic projection centered on Europe


# Start with the base map
map <- ggplot() +
  geom_sf(data = world_proj, fill = "#1C365F", color = "gray50") +
  # geom_sf_text(data = cities_proj, aes(label = name), nudge_y = 0) +
  coord_sf(crs = ortho_crs, 
           xlim = c(-6e6, 6e6), ylim = c(-3e6, 6e6), expand = FALSE) +
  theme_minimal()




tree_colors <- tibble(
  taxo=factor(c("TSAR","Haptista","Cryptista","Archaeplastida","Amorphea","Obazoa","Excavata","Eukaryota_X", "other")),
  color=c("#b25545","#b08f89", "#8c8fc6", "#afb96c", "#959ca5","#959ca5", "#917990ff", "#ac9e71ff", "grey")
)
w <- 7e5  # ~300 km
h <- 15e4


# Annotate plots
for (i in seq_len(nrow(cities_coords))) {
  x <- cities_coords[i, "X"]
  y <- cities_coords[i, "Y"]
  
  plot_grob <- ggplotGrob(
    p_list[[i]] + geom_col(show.legend = FALSE) + 
      scale_fill_manual(breaks = tree_colors$taxo, values=tree_colors$color)+
      theme_void()
  )
  
  map <- map + annotation_custom(
    grob = plot_grob,
    xmin = x - w/2,
    xmax = x + w/2,
    ymin = y - h/2,
    ymax = y + h/2
  )
}


legend_grob <- cowplot::get_legend(p_list[[1]]+geom_col()+
                                     scale_fill_manual(breaks = tree_colors$taxo, values=tree_colors$color)+theme(legend.title = element_blank()))


full_p <- map+ annotation_custom(
  grob = legend_grob,
  xmin = -5e6, xmax = -15,  # Choose a good spot
  ymin = 48, ymax = 60
)+
  geom_sf_text(data = cities_proj, aes(label = name), nudge_y = 0, color="white")






pdf(file="/g/schwab/marco/projects/trec_metaB/europe_ov_plot_newproj.pdf", width=30, height=30)
print(full_p)
dev.off()




png(file="/g/schwab/marco/projects/trec_metaB/europe_ov_plot_newproj.png", width=5000, height=5000, res = 200)
print(full_p)
dev.off()



