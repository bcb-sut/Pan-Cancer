library(readr)
library(dplyr)
library(tidyr)

files = list.files("data/codings/", full.names = FALSE)

for (file in files) {
  cancer_type = sub(".tsv.gz", "", file)
  print(cancer_type)
  data = read_tsv(paste("data/codings/", file, sep = "/"))
  
  data %>% 
    select(icgc_sample_id, chromosome, position, strand, gene_id, motif_3mer) %>% 
    unique() -> mutations
             
  
  mutations %>% 
    group_by(icgc_sample_id, gene_id) %>%
    summarise(count = n()) %>%
    spread(gene_id, count, fill = 0) -> genes_mat

  # data %>% 
  #   # filter(gene_biotype == "protein_coding") %>% 
  #   dplyr::select(icgc_sample_id, gene_id) %>% 
  #   unique() %>%
  #   mutate(count = 1) %>% 
  #   spread(gene_id, count, fill = 0) -> gene_mat
  
  mutations %>% 
    mutate(genemotif = paste(gene_id, motif_3mer, sep = "_")) %>% 
    group_by(icgc_sample_id, genemotif) %>% 
    summarise(count = n()) %>%
    spread(genemotif, count, fill = 0) -> genemotifs_mat
  
  write_tsv(genes_mat,
            file.path("matrices/allgenes/", paste0(cancer_type, ".tsv.gz")))
  write_tsv(genemotifs_mat,
            file.path("matrices/allgenemotifs/", paste0(cancer_type, ".tsv.gz")))
}




