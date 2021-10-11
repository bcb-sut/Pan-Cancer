library(readr)
library(dplyr)
library(tidyr)

gene_fantoms = list.files(path = "data/gene_fantom/", 
                   full.names = FALSE, recursive = FALSE)
motifs = list.files(path = "data/motifs/", 
                   full.names = FALSE, recursive = FALSE)

n = length(motifs)

for (i in 1:n) {
  cancer_type = sub(".tsv", "", gene_fantoms[i])
  print(cancer_type)
  
  gene_df = read_tsv(paste0("data/gene_fantom/", gene_fantoms[i]))
  motif_df = read_tsv(paste0("data/motifs/", motifs[i]))
  
  motif_df %>%
    dplyr::select(icgc_sample_id, chromosome, position, strand, motif_3mer, motif_5mer) %>% unique() -> motif_df

  gene_df %>%
    filter(gene_biotype == "protein_coding") %>%
    dplyr::select(icgc_sample_id, chromosome, position, strand, gene_id) %>%
    left_join(motif_df) %>% 
    drop_na() %>% 
    mutate(genemotif_3mer = paste(gene_id, motif_3mer, sep = "_")) %>% 
    mutate(genemotif_5mer = paste(gene_id, motif_5mer, sep = "_")) %>% 
    dplyr::select(-c(gene_id, motif_3mer, motif_5mer)) -> genemotif_df
  
  write_tsv(genemotif_df, file.path("data/genemotifs/", paste0(cancer_type, "_genemotif.tsv.gz")))

}

