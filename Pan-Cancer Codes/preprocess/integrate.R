library(readr)
library(tidyr)
library(dplyr)

hg19 = read_tsv("data/ENCODE_hg19.tsv")

hg19 %>% 
  rename(gene_id = gene_name) %>% 
  mutate(gene_id = sub("\\.\\d+", "", gene_id)) %>% 
  select(gene_id, Gene_class) -> hg19

files = list.files(path = "data/motifs/", 
                    full.names = FALSE, recursive = FALSE)

for (file in files) {
  cancer_type = sub("_motif.tsv", "", file)
  print(cancer_type)
  data = read_tsv(file.path("data/motifs", file))
  
  data %>% 
    rename(gene_id = gene_affected) %>% 
    # rename(sample_id = icgc_sample_id) %>%
    left_join(hg19, by = "gene_id") %>% 
    filter(Gene_class == "protein_coding") %>% 
    select(icgc_sample_id, chromosome, position, strand, gene_id, motif_3mer, motif_5mer, sequencing_strategy) -> final_data
  
  write_tsv(final_data, file.path("data/codings/", paste0(cancer_type, ".tsv.gz")))
  
}