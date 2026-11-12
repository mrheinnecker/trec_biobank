


library(tidyverse)
library(getopt)
library(googlesheets4)
library(googledrive)
library(cowplot)
#email = "marco.rheinnecker@embl.de"

json_key <- "/g/schwab/marco/repos/tem_classification/tem/trec-tem-screen-e98a2e03f58b.json"
gs4_auth(path=json_key)
drive_auth(path = json_key)
trec_tem_googledoc <- "https://docs.google.com/spreadsheets/d/143uVeeJ72SQE5eK01lzWYCEiT7pJUF3lX7hJl3R9s9I/edit?gid=258669282#gid=258669282"



df_trec_tem_current_state <- read_sheet(trec_tem_googledoc, sheet="image_log", col_types="c") 


dic <- tibble(
  old=c("freezing_damage" , "infiltration_issue", "knife_marks", "folds", "precipitation", "biological_artifacts", "imaging_artifacts","blending_artifacts" ),
  new=c("Freezing damage", "Infiltration artifacts", "Knife marks", "Folds", "Post-staining artifacts",
        "Ultrastructural integrity", "Imaging artifacts","Blending artifacts")
)


annotated <- df_trec_tem_current_state %>%
  filter(!is.na(user), dinoflagellate==0) 

all_sites <- c("TAL","KRI","BIL","POR","BAR","NAP","ATH") #%>% c(., "VIG")

df <- annotated %>%
  select(shortname, site,  blending_artifacts,imaging_artifacts,folds, freezing_damage, 
         precipitation, infiltration_issue, knife_marks, biological_artifacts) %>%
  pivot_longer(cols=names(.) %>% .[which(!. %in% c("shortname", "site"))]) %>%
  mutate(
    value=ifelse(is.na(value), 0, value),
    val_tf=ifelse(value==0|is.na(value), 0, 1),
    site=factor(site, levels = all_sites),
    name=factor(name, levels=rev(dic$old))
  )

cell_count <- annotated %>% mutate(site=factor(site, levels = all_sites)) %>% group_by(site, .drop=F) %>% tally() %>% select(site, total_cells=n)

plot_data <- 
  df %>%
  filter(val_tf!=0) %>% 
  group_by(site, name, .drop=F) %>%
  tally() %>%
  left_join(
    cell_count
  ) %>%
  mutate(
    frac=n/total_cells
  )
  
  
  

hm <- ggplot(plot_data, aes(x=site, fill=frac, y=name))+
  geom_tile()+
  geom_text(aes(label=n), color="white")+
  scale_fill_gradientn(colors=c("#FFC20A","#0C7BDC"), name="fraction")+
  scale_y_discrete(breaks=dic$old, labels=dic$new)+
  theme_bw()+
  ylab("artifact")



bp <- ggplot(cell_count, aes(x=site, y=total_cells))+
  geom_col()+
  geom_text(aes(label=total_cells, y=total_cells-1), color="white", vjust=1)+
  geom_label(data=tibble(x=ifelse("VIG" %in% all_sites, "ATH", "NAP"), y=130, label=paste0("total: ", as.character(sum(cell_count$total_cells)))),
             aes(x=x, y=y, label=label), inherit.aes=F)+
  theme_bw()+
  ylab("total cells")


comb_plot <- cowplot::plot_grid(
  bp+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()), 
  hm,
  align="v",
  rel_heights=c(1.5,3),
  nrow=2
)


outdir <- "/g/schwab/tem_screen"

plot_name <- "TEM_screen_artifact_count"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=5.5, height=4.5)
comb_plot
dev.off()




png(
  filename = file.path(outdir, paste0(plot_name, ".png")),
  width = 5.5,
  height = 4.5,
  units = "in",
  res = 300   # good default for publication-quality
)

comb_plot
dev.off()





remove_level <-  
  df %>%
  mutate(
    name=factor(name, levels=dic$old),
    value=factor(value)
  )%>%
  group_by(name, value, .drop=F
  ) %>%
    tally() %>%
    filter(n==0) %>%
  mutate(remove=T) %>%
  select(-n)



dpd <- df %>%
  mutate(
    name=factor(name, levels=dic$old),
    value=factor(value)
  )%>%
  group_by(site, name, value, .drop=F
           ) %>%
  tally() %>%
  #filter(!name %in% c("imaging_artifacts", "blending_artifacts")) %>%
  left_join(remove_level) %>%
  filter(is.na(remove))

#dpd %>%
#  filter(value==0, n==0)


## imaging and knife marks rasunehmen -> 6 panels
## color code match with main figure
## to do: write imagin processing paragraph



cols <- c("#FFC20A", "#A3A500", "#4DAF7C",  "#0C7BDC", "#6A3D9A")


lapply(unique(dpd$name), function(ART){
  
  print(ART)
  
  rel_dpd <- dpd %>% filter(name==ART)
    
  pd <- position_dodge(width = 0.9)
  
  detailed_plot <- ggplot(rel_dpd, aes(x = site, y = n, fill = value)) +
    #facet_wrap(~site, scales = "free_y", ncol = 2) +
    geom_col(
      position=pd
      #position = position_dodge2(width = 0.9, preserve = "single"),
      #width = 1
    ) +
    
    #scale_x_discrete(breaks=dic$old, labels=dic$new)+
    scale_fill_manual(name="level", values=rev(cols))+
    theme_bw() +
    #theme(axis.text.x = element_text(angle = 325, hjust = 0))+
    xlab("site")+
    ylab("number of occurrences")
  
  
  
  plot_name_ls <- paste0(ART, "_landscape")
  
  pdf(file=file.path(outdir, "detailed_subplots",paste0(plot_name_ls, ".pdf")), width=4, height=3)
  print(detailed_plot+geom_text(aes(label=n), position=pd, vjust=0, color="black", size=2.5))
  dev.off()
  
  png(
    filename = file.path(outdir, "detailed_subplots", paste0(plot_name_ls, ".png")),
    width = 4,
    height = 3,
    units = "in",
    res = 300   # good default for publication-quality
  )
  
  print(detailed_plot+geom_text(aes(label=n), position=pd, vjust=0, color="black", size=2.5))
  dev.off()
  
  
  plot_name_pt <- paste0(ART, "_portrait")
  
  pdf(file=file.path(outdir, "detailed_subplots",paste0(plot_name_pt, ".pdf")), width=4, height=3)
  print(detailed_plot+geom_text(aes(label=n), position=pd, hjust=0, color="black", size=2.5)+coord_flip())
  dev.off()
  
  png(
    filename = file.path(outdir, "detailed_subplots", paste0(plot_name_pt, ".png")),
    width = 4,
    height = 3,
    units = "in",
    res = 300   # good default for publication-quality
  )
  
  print(detailed_plot+geom_text(aes(label=n), position=pd, hjust=0, color="black", size=2.5)+coord_flip())
  dev.off()
  
  
  
})








