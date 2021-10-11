library(readr)
library(dplyr)
library(tidyr)

files = list.files(path = "data/codings/", 
                   full.names = FALSE, recursive = FALSE)
View(files)
genes_mat = data.frame(matrix(ncol = 3, nrow = 0))
colnames(genes_mat) <- c("gene_id", "count", "cancer_type")

for(file in files) {
  cancer_type = sub(".tsv.gz", "", file)
  print(cancer_type)
  data = read_tsv(paste0("data/codings/", file))
  View(data)
  data %>% 
    filter(!is.na(gene_id)) %>%
    # filter(gene_biotype == "protein_coding") %>%
    dplyr::select(icgc_sample_id, gene_id) %>% 
    unique() %>%
    group_by(gene_id) %>% 
    dplyr::summarise(count = n()) %>%
    mutate(cancer_type = cancer_type) -> patients_count
  # write_tsv(genes,
            # file.path("gene_count/", "gene_count.tsv"))
  genes_mat = rbind(genes_mat, patients_count)
}

write_tsv(genes_mat, file.path("data/counts/gene_count.tsv.gz"))
