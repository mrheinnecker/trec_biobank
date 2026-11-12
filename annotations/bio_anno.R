


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



annotated_bio <- df_trec_tem_current_state %>%
  filter(!is.na(`non-identified`)) %>%
  select(shortname, site,
         theca, 
         nucleus, 
         nucleolus, 
         golgi_apparatus, 
         mitochondria, 
         chloroplasts, 
         plastoglobuli, 
         starch, 
         vacuole, 
         ER, 
         tubular_network,
         retic_net_cryst_net=`reticulated_network/crystalline_compartment`, 
         electron_dense_sheets, 
         pusule, 
         extrusomes,  
         rhabdosome, 
         eyespot, 
         flagellum=`basal_body/flagellum`, 
         symbiosis_parasitism=`symbiotic/parasitism`,
         non_identified=`non-identified`)



colsorter <- names(annotated_bio)[3:22]



col_high <- "#0C7BDC"
col_low <- "#FFC20A"



pd <- annotated_bio %>%
  pivot_longer(cols=names(.) %>% .[which(!. %in% c("shortname", "site"))])%>%
  mutate(
         name=factor(name, levels=colsorter)) %>%
  filter(value!=0) %>%
  group_by(name, .drop=F) %>%
  tally() %>%
  arrange(n) %>%
  mutate(
    frac=n/nrow(annotated_bio)
  )



p_abundance <- ggplot(pd, aes(x=name, y=n))+
  geom_col(aes(fill=frac))+
  geom_text(aes(label=n), vjust=0)+  scale_fill_gradientn(
    colors = c("grey", col_low, col_high),
    values = scales::rescale(c(0, 0.01, 1)),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    name="fraction"
  )+
  theme_bw()+
  #facet_wrap(~site)+
  theme(axis.text.x = element_text(angle=325, hjust=0))+
  ylab("total cells")


outdir <- "/g/schwab/tem_screen"

plot_name <- "TEM_screen_organelle_abundance"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=7, height=4.5)
p_abundance
dev.off()

all_sites <- c("TAL","KRI","BIL","POR","BAR","NAP","ATH")

pdhm <- annotated_bio %>%
  pivot_longer(cols=names(.) %>% .[which(!. %in% c("shortname", "site"))]) %>%
  filter(value!=0) %>%
  mutate(site=factor(site, levels = all_sites), name=factor(name, levels=colsorter)) %>%
  group_by(site, name, .drop=F) %>%
  tally() %>%
  left_join(annotated_bio %>% group_by(site) %>% tally() %>% select(n_per_site=n, site)) %>%
  mutate(frac=n/n_per_site)





p_hm <- ggplot(pdhm, aes(x=name, y=site, fill=frac))+
  geom_tile()+
  geom_text(aes(label=n), color="white")+
  
  scale_fill_gradientn(
    colors = c("grey", col_low, col_high),
    values = scales::rescale(c(0, 0.01, 1)),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    name="fraction"
  )+
  theme_bw()+
  #facet_wrap(~site)+
  theme(axis.text.x = element_text(angle=325, hjust=0))+
  xlab("cell structure / organelle")


outdir <- "/g/schwab/tem_screen"

plot_name <- "TEM_screen_organelle_abundance_per_site"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=7, height=4.5)
p_hm
dev.off()



comb_plot <- cowplot::plot_grid(
  p_abundance+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "none"), 
  p_hm,
  align="v",
  rel_heights=c(2,3),
  nrow=2
)


plot_name <- "TEM_screen_organelle_abundance_combined"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=10, height=7.5)
comb_plot
dev.off()





png(
  filename = file.path(outdir, paste0(plot_name, ".png")),
  width = 10,
  height = 7.5,
  units = "in",
  res = 300   # good default for publication-quality
)

comb_plot
dev.off()




## fisher tests


cloro_starch <- fisher.test(annotated_bio$chloroplasts, 
            annotated_bio$starch)


cloro_tub <- fisher.test(annotated_bio$tubular_network, 
                            annotated_bio$chloroplasts)


dct <- annotated_bio %>%
  group_by(chloroplasts, tubular_network) %>%
  tally() 


dcs <- annotated_bio %>%
  group_by(chloroplasts, starch) %>%
  tally() 


pct <- ggplot(dct, aes(y=chloroplasts, x=tubular_network, fill=n))+
  geom_tile(show.legend=F)+
  geom_text(aes(label=n), color="white")+
  geom_text(inherit.aes = F, color="black", vjust=0,
            data=tibble(or=cloro_tub$estimate, p=cloro_tub$p.value), 
            aes(x=1.5, y=2.6, label=paste("OR:", round((or),2), "\np:", round(p,5))))+
  scale_fill_gradientn(colors = rev(c(col_high, "grey80")))+
  scale_x_discrete(breaks=c(0,1), labels=c("False", "True"), expand = expansion(mult = c(0.5, 0.5)))+
  scale_y_discrete(breaks=c(0,1), labels=c("False", "True"), expand = expansion(mult = c(0.5, 1.3)))+
  #scale_y_
  theme_bw()+
  xlab("tubular network")+
  theme(panel.grid = element_blank())

pcs <- ggplot(dcs, aes(y=chloroplasts, x=starch, fill=n))+
  geom_tile(show.legend=F)+
  geom_text(aes(label=n), color="white")+
  geom_text(inherit.aes = F, color="black", vjust=0,
            data=tibble(or=cloro_starch$estimate, p=cloro_starch$p.value), 
            aes(x=1.5, y=2.6, label=paste("OR:", round((or),2), "\np:", round(p,5))))+
  scale_fill_gradientn(colors = rev(c(col_high, "grey80")))+
  scale_x_discrete(breaks=c(0,1), labels=c("False", "True"), expand = expansion(mult = c(0.5, 0.5)))+
  scale_y_discrete(breaks=c(0,1), labels=c("False", "True"), expand = expansion(mult = c(0.5, 1.3)))+
  theme_bw()+
  theme(panel.grid = element_blank())


comb_plot <- cowplot::plot_grid(
  pct, pcs
)





plot_name <- "fisher_tests"

pdf(file=file.path(outdir, paste0(plot_name, ".pdf")), width=5, height=2)
comb_plot
dev.off()



png(
  filename = file.path(outdir, paste0(plot_name, ".png")),
  width = 5,
  height = 2,
  units = "in",
  res = 300   # good default for publication-quality
)

comb_plot
dev.off()





















