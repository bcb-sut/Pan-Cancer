library(readr)
library(dplyr)

allgenes = read_tsv("matrices/allgenes.tsv.gz")

sub_mat <- function(sgnf_genes) {
  sgnf_genes = sgnf_genes$gene_id
  names.use = names(allgenes)[names(allgenes) %in% sgnf_genes]
  names.use = c("icgc_sample_id", names.use)
  mat = allgenes[, names.use]
  return(mat)
}

mat_path = "matrices/sample_sgnfgene/"

sgnf_genes = read_tsv("sgnf_genes/allcancer_0.005.tsv")
write_tsv(sub_mat(sgnf_genes), file.path("matrices/sgnf_genes/", "allcancer_0.005.tsv.gz"))

