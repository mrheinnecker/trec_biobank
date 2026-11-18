
#query <- "Heterocapsa"
#query <- "Dinoflag"


europe_ov_plot <- function(query){
  
  
  
}

taxonomy_mapping <- function(){
  
  return(
    tibble(
      level=paste0("level_", 1:9),
      pr2_taxonomy_levels = c(
        "Domain",
        "Supergroup",
        "Division",
        "Subdivision",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      )
      
    )
  )
  
}


get_highest_score_per_level <- function(query, df_asv_taxonomy, highest=1){
  
  level_cols <- df_asv_taxonomy %>% select(starts_with("level_"))
  
  hit_per_level <- lapply(names(level_cols), function(COLNAME){
    
    #print(level_cols[[COLNAME]])
    
    full_list <- level_cols[[COLNAME]] %>% unique() %>% as.character() 
    
    distances <- tibble(entry=full_list,
                        score=stringdist::stringdist(query, full_list, method = "jw"))
    
    hit <- distances %>%
        arrange(score) %>%
        .[1:highest,] %>%
      mutate(level=COLNAME)
    
    return(hit)
  }) %>%
    bind_rows()
  
  return(hit_per_level)

}




get_highest_score_colname <- function(query, df_asv_taxonomy){
  
  hit_per_level <- get_highest_score_per_level(query, df_asv_taxonomy)
  rel_col <- which(hit_per_level$score==min(hit_per_level$score))
  colname <- paste0("level_", rel_col)
  
}

tree_taxo_new <- function(query, df_asv_taxonomy){
  hit_per_level <- get_highest_score_per_level(query, df_asv_taxonomy, 3)


  p_text_hits <- ggplot(hit_per_level)+
    geom_text(
      aes(x=level, y=-score, label=entry, key=entry)
    )+
    scale_x_discrete(breaks=taxonomy_mapping() %>% pull(level),
                     labels=taxonomy_mapping() %>% pull(pr2_taxonomy_levels))+
    
    theme_bw()+
    theme(
      axis.text.x = element_text(angle=315),
      axis.text.y = element_blank()
    )+
    ylab("similarity")

  return(p_text_hits)
}


prepare_tree_taxo_data <- function(query, df_asv_taxonomy){
  
  hit_per_level <- get_highest_score_per_level(query, df_asv_taxonomy)
  
  rel <- hit_per_level[which(hit_per_level$score==min(hit_per_level$score)),]
  
  rel_col <- rel$level %>% str_split("_") %>% map_chr(.,2) %>% as.numeric()
  
  colname <- rel$level
  hit <- rel$entry
  
  
  sel_cols <- paste0("level_", seq(1, rel_col+1))
  
  plot_data_raw <- 
    df_asv_taxonomy %>%
    filter(select(., all_of(colname))==hit) %>%
    select(all_of(sel_cols)) %>%
    unique()%>%
    pivot_longer(cols=names(.), names_to = "level", values_to = "taxo") %>%
    mutate(num_level=str_replace(level, "level_", "") %>% as.numeric()) %>%
    unique() %>%
    group_by(num_level) %>%
    mutate(ypos_raw=seq(1,length(num_level))) 
  
  return(plot_data_raw)
  
}


tree_taxo <- function(query, df_asv_taxonomy){
  
  
  # target_list <- df_asv_taxonomy %>% select(starts_with("level_")) %>%
  #   as.character() %>% unique()
  # 
  # 
 # test_mat <- matrix(target_list, nrow=4)
  
  # hit_per_level <- get_highest_score_per_level(query, df_asv_taxonomy)
  # 
  # rel <- hit_per_level[which(hit_per_level$score==min(hit_per_level$score)),]
  # 
  # rel_col <- rel$level %>% str_split("_") %>% map_chr(.,2) %>% as.numeric()
  # 
  # colname <- rel$level
  # hit <- rel$entry
  # 
  # 
  # sel_cols <- paste0("level_", seq(1, rel_col+1))
  # 
  # plot_data_raw <- 
  #   df_asv_taxonomy %>%
  #   filter(select(., all_of(colname))==hit) %>%
  #   select(all_of(sel_cols)) %>%
  #   unique()%>%
  #   pivot_longer(cols=names(.), names_to = "level", values_to = "taxo") %>%
  #   mutate(num_level=str_replace(level, "level_", "") %>% as.numeric()) %>%
  #   unique() %>%
  #   group_by(num_level) %>%
  #   mutate(ypos_raw=seq(1,length(num_level))) 
  
  plot_data_raw <- prepare_tree_taxo_data(query, df_asv_taxonomy)
  #print(plot_data_raw)
  
  plot_data <- plot_data_raw %>%
    left_join(
      plot_data_raw %>% group_by(num_level) %>% tally()
    ) %>%
    ungroup() %>%
    mutate(
      ypos=ifelse(n==1, 
                  0.5*(max(plot_data_raw$ypos_raw)+1),
                  ypos_raw),
      ang=ifelse(n==1, 
                 315,
                 0),
      just=ifelse(n==1, 
                  1,
                  0),
      lab=ifelse(n==1, 
                 paste0(taxo, " "),
                 paste0(" ", taxo)),
    )
  
  line_data <- plot_data %>%
    filter(num_level %in% c(max(plot_data$num_level))) %>%
    select(y=ypos, x=num_level) %>%
    mutate(yend=plot_data %>%
             filter(num_level %in% c(max(plot_data$num_level)-1)) %>% pull(ypos),
           xend=max(plot_data$num_level)-1) %>%
    bind_rows(
      tibble(x=1, xend=max(plot_data$num_level)-1,
             y=0.5*(max(plot_data_raw$ypos_raw)+1),
             yend=0.5*(max(plot_data_raw$ypos_raw)+1))
    )
  
  tree_plot <- ggplot(plot_data,
         aes(x=num_level, y=ypos))+
    geom_segment(data=line_data, aes(x=x, y=y, xend=xend, yend=yend))+
    geom_point(color="red", aes(key=taxo))+
    geom_text(aes(label=lab, angle=ang, hjust=just))+
    theme_void()
  
  return(tree_plot)
  
}





europe_map_plot <- function(){
  
  
  
  
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
  
  
  
  return(full_p)
  
}


sequence_query <- function(query_seq, df_asv_taxonomy, similarity=0.9){
  
  
  all_seqs <- DNAStringSet(df_asv_taxonomy$sequence, use.names = T)
  names(all_seqs) <- df_asv_taxonomy$asv_id
  
  #query_seq <- df_asv_taxonomy$sequence[4000] 
  
  query_seq <-"GGTGAAAGCCCATCGCTCAACGGTGGAACGGCCATTGATACTGTCTGACTTGAATTATTAGGAAGTAACTAGAATATGTAGTGTAGCGGTGAAATGCTTAGAGATTACATGGAATACCAATTGCGAAGGCAGGTTACTACTAATTGATTGACGCTGATGGACGAAAGCGTGGGTAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGATACTAGCTGTTGGGGGCAACTTCAGTGGCTAAGCGAAAGTGATAAGTATCCCACCTGGGGAGTACGTTCGCAAGAATG"
  
  db_matrix <- oligonucleotideFrequency(DNAStringSet(all_seqs), width=5)
  query_vector <- oligonucleotideFrequency(DNAString(query_seq), width=5)
  db_matrix <- as.matrix(db_matrix)
  query_vector <- as.matrix(query_vector)
  
 # rownames(db_matrix) <- paste0("seq_", seq_len(nrow(db_matrix)))
  ## maybe this can be accelerated
  sim_scores <- apply(db_matrix, 1, function(row_vec){
    
    cosine(as.numeric(query_vector), row_vec)
    
  })
  
  
  max(sim_scores)
  
  names(all_seqs[which(sim_scores>0.5)])
  
}


fill_up_full_hm_grid <- function(asvs_per_sample_frac){
  
  sample_site_map <- asvs_per_sample_frac %>%
    select(sample_id, site) %>%
    distinct()
  
  # Then, create the full grid
  full_grid <- expand.grid(
    sample_id = unique(sample_site_map$sample_id),
    asv_id = unique(asvs_per_sample_frac$asv_id),
    stringsAsFactors = FALSE
  ) %>%
    left_join(sample_site_map, by = "sample_id")
  
  # Then join the original data and fill NAs with 0
  plot_data <- full_grid %>%
    left_join(
      asvs_per_sample_frac %>% 
        filter(asv_id %in% asvs_per_sample_frac$asv_id) %>%
        select(sample_id, asv_id, frac),
      by = c("sample_id", "asv_id")
    ) %>%
    mutate(frac = replace_na(frac, 0)) %>%
    as_tibble()
  
  return(plot_data)
  
}


select_taxo_level <- function(level_inputs, df_asv_per_sample, selected_sample=NULL){
  
  first_non_selected_level <- min(c(which(level_inputs=="--all--"),9))
  print(paste("levels:", first_non_selected_level))
  # Dummy plot showing taxa counts (replace with real logic if needed)
  if(is.null(selected_sample)){
    df_filtered <- df_asv_per_sample %>%
      left_join(df_asv_taxonomy) 
  } else {
      df_filtered <- df_asv_per_sample %>%
    filter(sample_id == selected_sample,
    ) %>%
    left_join(df_asv_taxonomy) 
  }

  
  #print(names(df_filtered))
  
  for (LEV in seq(2,first_non_selected_level)){
    
    df_filtered <- df_filtered[which(df_filtered[[paste0("level_", LEV-1)]] == level_inputs[[LEV-1]]),] 
    
    #print(nrow(df_filtered))
    
  }
  
  df_filtered$col_x <- df_filtered[[paste0("level_", first_non_selected_level)]]
  df_filtered$col_fill <- df_filtered[[paste0("level_", first_non_selected_level+1)]]

  
  return(df_filtered)
  
}


get_main_abundance_plot <- function(asvs_per_sample, df_asv_taxonomy, df_asv_per_sample, 
                                    norm_type="Eukaryota", ylab_in="asv_id", grouping="asv_id"){
  mx_frac <- 1
  legend_scale <- c(0,0.01,0.5,round(mx_frac*20)/20)
  
  
  color_vec <- c("black","blue4","blue3","darkcyan","green4","green","yellow4","yellow","pink", "deeppink")
  val_vec <- c(0,
               0.001/mx_frac,
               0.0099/mx_frac,
               0.01/mx_frac,
               0.0999/mx_frac,
               0.25/mx_frac,
               0.4999/mx_frac,
               0.5/mx_frac,
               0.74999/mx_frac,
               0.75/mx_frac,
               1)
  
  asvs_to_norm <- df_asv_taxonomy %>%
    filter(level_1==norm_type) %>%
    pull(asv_id)
  
  tot_sample <- df_asv_per_sample %>%
    filter(asv_id %in% asvs_to_norm)%>%
    group_by(sample_id) %>%
    summarize(
      to_norm_reads=sum(nreads)
    )
  
  print(1)
  summed_barplot_data <- asvs_per_sample %>%
    group_by(sample_id) %>%
    summarize(
      tot_reads=sum(nreads),
      site=unique(site)
    ) %>%
    left_join(
      tot_sample, by="sample_id"
    ) %>%
    mutate(
      frac=tot_reads/to_norm_reads
    )
  
  
   export_priorities <- summed_barplot_data  %>%
     arrange(site, desc(frac)) %>%
     group_by(site) %>%
     mutate(
       priority=c(1:length(site))
     ) %>%
     filter(
       priority<=500
     ) %>%
     left_join(df_sampling_sites %>% select(sample_id, date, time, size_fraction))
     
   write_tsv(export_priorities, file="/home/rheinnec/export_priorities.tsv")
  # 
  # 
  print(2)
  bar_sumplot <- ggplot(summed_barplot_data,
         aes(x=as.character(sample_id), y=tot_reads, fill=frac))+
    facet_grid(~site, scales = "free_x", space = "free")+
    geom_col()+
    scale_fill_gradientn(
      # name="abundance [%]",
      # values=val_vec,
      # breaks=legend_scale,
      #labels=paste(legend_scale*100),
      colors = color_vec
    )+
    theme_bw()+
    theme(
          axis.text.x = element_text(angle = 270, hjust=0.5, vjust=0.5),
        
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  print(3)
  if(grouping=="asv_id"){
    pregrouped <- asvs_per_sample
  } else if(grouping=="level_9"){
    pregrouped <- asvs_per_sample %>%
      select(asv_id, sample_id, nreads, site) %>%
      left_join(df_asv_taxonomy) %>%
      group_by(sample_id, level_9) %>%
      summarize(
        nreads=sum(nreads),
        site=unique(site),
        
        #sample_id=sample_id
      ) %>%
      ungroup() %>%
      mutate(asv_id=level_9)
  } else if(grouping=="level_8"){
    pregrouped <- asvs_per_sample %>%
      select(asv_id, sample_id, nreads, site) %>%
      left_join(df_asv_taxonomy) %>%
      group_by(sample_id,level_8) %>%
      summarize(
        nreads=sum(nreads),
        site=unique(site),
        
        #sample_id=unique(sample_id)
      )%>%
      ungroup()%>%
      mutate(asv_id=level_8)
  } else if(grouping=="level_7"){
    pregrouped <- asvs_per_sample %>%
      select(asv_id, sample_id, nreads, site) %>%
      left_join(df_asv_taxonomy) %>%
      group_by(sample_id,level_7) %>%
      summarize(
        nreads=sum(nreads),
        site=unique(site),
        
        #sample_id=unique(sample_id)
      )%>%
      ungroup()%>%
      mutate(asv_id=level_7)
  } else if(grouping=="level_6"){
    pregrouped <- asvs_per_sample %>%
      select(asv_id, sample_id, nreads, site) %>%
      left_join(df_asv_taxonomy) %>%
      group_by(sample_id,level_6) %>%
      summarize(
        nreads=sum(nreads),
        site=unique(site),
        
        #sample_id=unique(sample_id)
      )%>%
      ungroup()%>%
      mutate(asv_id=level_6)
  }
  
  
  asvs_per_sample_frac <- pregrouped %>%
    group_by(sample_id) %>%
    mutate(tot_sample_reads=sum(nreads)) %>%
    rowwise() %>%
    mutate(
      frac=nreads/tot_sample_reads
    )
  
  print(4)
  plot_data <- fill_up_full_hm_grid(asvs_per_sample_frac)
  
  print(5)
  ylabs <- plot_data %>%
    left_join(df_asv_taxonomy, by="asv_id")


  ylabs$ylab <- ylabs[[ylab_in]]
  
  print(6)
  date <- Sys.Date()
  time <- Sys.time()
  write_tsv(plot_data, file=paste0("~/plot_data", date,time,  ".tsv"))
  
  hm <- 
    ggplot(plot_data# %>%
               #  mutate(asv_id = factor(asv_id, levels = asv_phylo_order))
               )+
    geom_tile(
      aes(x=as.character(sample_id), y=asv_id, fill=frac),
      #show.legend = F
    )+
    facet_grid(~site, scales = "free_x", space = "free")+
    scale_fill_gradientn(
      name="abundance [%]",
      values=val_vec,
      breaks=legend_scale,
      labels=paste(legend_scale*100),
      colors = color_vec
    )+
    scale_y_discrete(breaks=ylabs$asv_id, labels=ylabs$ylab)+
    # scale_y_discrete(name=paste0("ASV [", annotations[which(annotations$name==y_annotation),]$labels, "]"),
    #                  breaks=dino_high_abundance$asv_id, labels=dino_high_abundance[[y_annotation]])+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 270, hjust=0.5, vjust=0.5)
    )
    # theme(
    #   legend.position = "bottom",
    #   axis.text = element_blank(),
    #   axis.ticks.x = element_blank(),
    #   axis.ticks.y = element_blank(),
    #   axis.title.y = element_blank()
    # )
  
  
  combined_plot <- cowplot::plot_grid(
    bar_sumplot,
    hm,
    ncol=1,
    align = "v",
    axis = "lr",
    rel_heights=c(1,max(c(2, 0.05*length(unique(pregrouped$asv_id)))))
  )
  
  print("plot fully created... returning to main script")
  
  return(list(combined_plot, nrows=max(c(length(unique(pregrouped$asv_id)), 30))))
  
}




asv_overview_per_sample <- function(asv_frac_per_sample, dino_high_abundance, tree, y_annotation="level_7", tree_xlim=0.2){
  
  
  annotations <- tibble(
    
    level=c(0:8),
    name=paste("level", c(9:1), sep="_"),
    colors=c("purple4","purple", "deeppink", "red", "orange", "orange4", "brown","brown4", "black"),
    labels=c("species", "genus", "family", "order","class", "subdivision", "division", "supergroup", "domain")
    
  )
  
  
  rel_asvs_full <- asv_frac_per_sample %>%
    filter(asv_id %in% dino_high_abundance$asv_id) #%>%
    #dplyr::rename(sample_id=sample)
  # First, get the list of all relevant samples and their sites
  # sample_site_map <- rel_asvs_full %>%
  #   select(sample, site) %>%
  #   distinct()
  # 
  # # Then, create the full grid
  # full_grid <- expand.grid(
  #   sample = unique(sample_site_map$sample),
  #   asv_id = unique(dino_high_abundance$asv_id),
  #   stringsAsFactors = FALSE
  # ) %>%
  #   left_join(sample_site_map, by = "sample")
  # 
  # # Then join the original data and fill NAs with 0
  # plot_data <- full_grid %>%
  #   left_join(
  #     rel_asvs_full %>% 
  #       filter(asv_id %in% dino_high_abundance$asv_id) %>%
  #       select(sample, asv_id, frac, level_4),
  #     by = c("sample", "asv_id")
  #   ) %>%
  #   mutate(frac = replace_na(frac, 0)) %>%
  #   as_tibble()
  plot_data <- fill_up_full_hm_grid(asv_frac_per_sample)
  
  
  mx_frac <- max(plot_data$frac)
  
  legend_scale <- c(0,0.01,0.05,0.1,0.2,round(mx_frac*20)/20)
   
  
  if(mx_frac>0.2){
    color_vec <- c("black","blue4","blue3","darkcyan","green4","green","yellow4","yellow","pink", "deeppink")
    val_vec <- c(0,
                 0.001/mx_frac,
                 0.0099/mx_frac,
                 0.01/mx_frac,
                 0.0499/mx_frac,
                 0.05/mx_frac,
                 0.0999/mx_frac,
                 0.1/mx_frac,
                 0.1999/mx_frac,
                 0.2/mx_frac,
                 1)
  } else {
    color_vec <- c("black","blue4","blue3","darkcyan","green4","green","yellow4","yellow")
    val_vec <- c(0,
                 0.001/mx_frac,
                 0.0099/mx_frac,
                 0.01/mx_frac,
                 0.0499/mx_frac,
                 0.05/mx_frac,
                 0.0999/mx_frac,
                 0.1/mx_frac,
                 1)
  }
  
  print(1)
  site_count <- rel_asvs_full %>% group_by(asv_id) %>% 
    summarise(all_sites=paste(unique(site), collapse = ";"),
              n_sites=length(unique(site)))
  
  asv_sorter <- plot_data %>%
    group_by(asv_id) %>%
    summarize(
      #n_sites=length(unique(site[which(frac>0)])),
      mn=mean(frac),
      md=median(frac),
      sm=sum(frac)
    ) %>%
    left_join(site_count) %>%
    arrange(desc(md))
  
  ## tree plot
  
  #library(patchwork)
  
  # Create tree plot (flip if you want to match heatmap)
  p_tree_raw <- 
    ggtree(tree) 
    #scale_x_continuous(limits=c(0,1.5))
  
  
  tree_data <- p_tree_raw$data
  
  # Get the tips in the order they appear (top to bottom by y position)
  asv_phylo_order <- tree_data %>%
    filter(isTip) %>%
    arrange(desc(y)) %>%  # use `desc(y)` if you want top-down order like in heatmaps
    pull(label)
  
  p_tree <- p_tree_raw + 
    scale_y_discrete(limits = asv_phylo_order)+
      coord_cartesian(xlim = c(0, tree_xlim))
  
  
 
  hm <- ggplot(plot_data %>%
                 mutate(asv_id = factor(asv_id, levels = asv_phylo_order)))+
    geom_tile(
      aes(x=as.character(sample_id), y=asv_id, fill=as.numeric(frac)),
      #show.legend = F
    )+
    facet_grid(~site, scales = "free_x", space = "free")+
    scale_fill_gradientn(
      name="abundance [%]",
      values=val_vec,
      breaks=legend_scale,
      labels=paste(legend_scale*100),
      colors = color_vec
    )+
    scale_y_discrete(name=paste0("ASV [", annotations[which(annotations$name==y_annotation),]$labels, "]"),
                     breaks=dino_high_abundance$asv_id, labels=dino_high_abundance[[y_annotation]])+
    theme_bw()+
    theme(
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  
  
  annotation_data <- dino_high_abundance %>%
    rowwise() %>%
    mutate(n_na_levels = sum(is.na(c_across(starts_with("level_"))))) %>%
    ungroup()
  
  
  site_count_plot <- 
    ggplot(asv_sorter %>% 
             left_join(annotation_data %>% select(asv_id, n_na_levels) %>% unique()) %>%
             mutate(asv_id=factor(asv_id, levels = asv_phylo_order)), 
           aes(y=asv_id, x=n_sites, fill=as.character(n_na_levels))) +
    geom_bar(stat="identity", #fill="lightblue"
             show.legend = T
    ) +
    #facet_grid(level_4~, scales = "free", space = "free")+
    # facet_wrap(~1)+
    #coord_flip() +
    theme_bw() +
    scale_fill_manual(
      name="annotation",
      breaks=annotations$level,
      values=annotations$colors,
      labels=annotations$labels)+
    scale_y_discrete(breaks=dino_high_abundance$asv_id, labels=dino_high_abundance[[y_annotation]],
                     position = "right")+
    theme(
      legend.position="bottom",
      #axis.text.y.right = element_text(),
      #axis.text.y.left = element_blank(),
      # axis.ticks.y = element_blank(),
      # axis.title.y = element_blank()
    )
  

  
  # Now use cowplot::plot_grid to combine the heatmap and the site count plot
  combined_plot <- 
    plot_grid(
      rel_widths = c(0.1, 1,0.25),
      p_tree,
      hm,           # Your original heatmap plot
      site_count_plot,  # The plot with total number of sites per ASV
      align = "h",  # Align plots horizontally
      ncol = 3,      # Arrange them in 2 columns
      axis = "bt"
    )
  
  return(combined_plot)
  
}


construct_phylogeny <- function(dino_high_abundance, df_asvs, method="nj"){
  
  
  
  all_cons_seqs <- dino_high_abundance %>%
    select(asv_id, level_6) %>%
    unique() %>%
    left_join(df_asvs %>% select(asv_id, sequence))
  
  
  
  dna_set <- DNAStringSet(all_cons_seqs$sequence)
  names(dna_set) <- all_cons_seqs$asv_id
  
  
  
  aligned_seqs <- AlignSeqs(dna_set)
  
  aln_mat <- as.matrix(aligned_seqs)
  
  
  # Identify columns where the proportion of gaps is below a threshold (e.g., keep positions with <50% gaps)
  keep_cols <- apply(aln_mat, 2, function(col) {
    mean(col == "-") < 0.5
  })
  
  # Subset the alignment
  trimmed_mat <- aln_mat[, keep_cols]
  
  # Convert back to DNAStringSet
  trimmed_seqs <- DNAStringSet(apply(trimmed_mat, 1, paste0, collapse=""))
  names(trimmed_seqs) <- names(aligned_seqs)
  
  if(method=="nj"){
    dm <- DistanceMatrix(trimmed_seqs, includeTerminalGaps=FALSE)
    tree <- nj(dm)
  } else {
    phy <- phyDat(as.matrix(trimmed_seqs), type = "DNA")
    
    # Create starting tree
    dm <- dist.ml(phy)
    treeNJ <- NJ(dm)
    
    # Fit a maximum likelihood model
    fit <- pml(treeNJ, data = phy)
    fit_opt <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic")
    
    tree <- fit_opt$tree
    
  }
  return(tree)
}



load_required_tables <- function(){
  
  mapping <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/v2/trec-aml-ssuv4v5_dada2_counts.tsv.gz")
  
  df_asvs <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/v2/trec-aml-ssuv4v5_dada2_asvs.tsv.gz")
  
  trec_id_mapping <- read_csv("/g/schwab/Chandni/MetaB_AML_40um/TREC_Plankton_MetaB_TRECAML_40um_dataset_files_export.csv") %>%
    mutate(
      sample=str_split(Name, "_") %>% map_chr(.,1),
      trec_name=str_replace_all(Samples, "SAMEA", "TREC")
    ) %>%
    select(err_id=sample, trec_id=Samples, trec_name) %>%
    unique()
  
  
  input_samples <- read_tsv("/g/schwab/marco/projects/trec_metaB/trec_metaB_site_table - Sheet1.tsv")%>%
    mutate(trec_id=paste0("SAMEA", Barcode_ID))
  
  
  error_samples <- read_xlsx("/g/schwab/marco/projects/trec_metaB/ErroneousReadset_ToRemove.xlsx")
  
  site_order <- c("Villefranche","Reykjavik","Roscoff", "Tallinn","Kristineberg", "Bilbao","Porto")
  
  site_mapping <- trec_id_mapping %>%
    left_join(input_samples) %>%
    select(err_id, site=Site) %>%
    mutate(
      site=case_when(
        str_detect(site, "Iceland") ~ "Reykjavik",
        str_detect(site, "Villef") ~ "Villefranche",
        TRUE ~ site
      ) %>% factor(., levels=site_order),
    )
  
  max_cols <- max(stringr::str_count(df_asvs$pr2_dada2_taxonomy, ";"), na.rm=T) + 1
  return(
    list(
      mapping=mapping,
      df_asvs=df_asvs %>%
        separate_wider_delim(pr2_dada2_taxonomy, delim = ";", 
                             names = paste0("level_", 1:max_cols), 
                             too_few = "align_start"),
      input_samples=input_samples,
      site_mapping=site_mapping,
      trec_id_mapping=trec_id_mapping,
      error_samples=error_samples
    )
  )
  
}

get_main_dataset <- function(df_asvs, site_mapping, mapping){
  
 
  
  taxo_anno <- df_asvs 
  
  
  
  
  plot_data_wide <- site_mapping %>% 
   
    #filter(site==SOI) %>%
    left_join(mapping, by=c("err_id"="sample")) %>%
    group_by(site, asv_id) %>%
    summarize(
      tot_reads=sum(nreads)
    ) %>%
    arrange(
      desc(tot_reads)
    ) %>%
    left_join(
      taxo_anno
    ) 
  
  return(plot_data_wide)
  
}
 


prepare_stacked_barplot_data <- function(df_asvs, site_mapping, mapping){
  
  
  
  plot_data <- get_main_dataset(df_asvs, site_mapping, mapping) %>%
    filter(level_1=="Eukaryota") %>%
    pivot_longer(cols = c("level_2"),
                 names_to = "level", values_to = "taxo") %>%
    
    
    
    ungroup()
  
  
  return(plot_data)
  
}

# df_asvs <- input_tables$df_asvs
# site_mapping <- input_tables$site_mapping
# mapping <- input_tables$mapping

get_stacked_barplot_per_site <- function(df_asvs, site_mapping, mapping, SITE){
  
  
 plot_data_all_sites <- prepare_stacked_barplot_data(df_asvs, site_mapping, mapping) %>%
   mutate(taxo=factor(taxo))
 
 
 taxo_levels <- levels(plot_data_all_sites$taxo)
 
 print(taxo_levels)
 
 plot_data <- plot_data_all_sites %>%
   filter(site==SITE)
 
 
 p <- ggplot(
   data=plot_data %>% group_by(taxo, site, .drop = F) %>% summarize(tot_reads_new=sum(tot_reads)) %>% 
     mutate(taxo=factor(taxo, levels=taxo_levels)) %>%
     filter(site==SITE),
   aes(
     y=site, x=tot_reads_new, fill=factor(taxo, levels = taxo_levels)
   )
 )
   #facet_wrap(~level, nrow=3)+
 
 
  
  return(p)
}










