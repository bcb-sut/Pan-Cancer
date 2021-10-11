library(readr)
library(dplyr)
library(tidyr)

files = list.files(path = 
                     "data/motifs/", 
                   full.names = FALSE, recursive = FALSE)

for(file in files) {
  data = read_tsv(paste0("data/motifs/", file))
  cancer_type = sub("_motif.tsv", "", file)
  print(cancer_type)
  data %>% 
    group_by(icgc_sample_id, motif_5mer) %>% 
    summarise(count = n()) %>% 
    spread(motif_5mer, count, fill = 0) %>% 
    ungroup() -> motif_mat_all
  data %>%
    filter(sequencing_strategy == "WGS") %>% 
    group_by(icgc_sample_id, motif_5mer) %>% 
    summarise(count = n()) %>% 
    spread(motif_5mer, count, fill = 0) %>% 
    ungroup() -> motif_mat_wgs
  data %>%
    filter(sequencing_strategy == "WXS") %>% 
    group_by(icgc_sample_id, motif_5mer) %>% 
    summarise(count = n()) %>% 
    spread(motif_5mer, count, fill = 0) %>% 
    ungroup() -> motif_mat_wxs
  
  write_path = "matrices/5mermotif/"
  
  write_tsv(motif_mat_wgs,
            file.path(write_path, paste0(cancer_type, "_WGS_motif.tsv")))
  
  write_tsv(motif_mat_wxs,
            file.path(write_path, paste0(cancer_type, "_WXS_motif.tsv")))
  
  write_tsv(motif_mat_all,
            file.path(write_path, paste0(cancer_type, "_motif.tsv")))
}


