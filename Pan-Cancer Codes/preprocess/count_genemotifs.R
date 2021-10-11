library(readr)
library(dplyr)
library(tidyr)

files = list.files(path = "data/codings/", 
                   full.names = FALSE, recursive = FALSE)

genes_mat = data.frame(matrix(ncol = 3, nrow = 0))
colnames(genes_mat) <- c("genemotif", "count", "cancer_type")

for(file in files) {
  cancer_type = sub(".tsv.gz", "", file)
  print(cancer_type)
  data = read_tsv(paste0("data/codings/", file))
  data %>% 
    filter(!is.na(gene_id)) %>%
    # filter(gene_biotype == "protein_coding") %>%
    mutate(genemotif = paste(gene_id, motif_3mer, sep = "_")) %>% 
    dplyr::select(icgc_sample_id, genemotif) %>% 
    unique() %>%
    group_by(genemotif) %>% 
    dplyr::summarise(count = n()) %>%
    mutate(cancer_type = cancer_type) -> patients_count
  # write_tsv(genes,
  # file.path("gene_count/", "gene_count.tsv"))
  genes_mat = rbind(genes_mat, patients_count)
}

write_tsv(genes_mat, file.path("data/counts/genemotif_count.tsv.gz"))

