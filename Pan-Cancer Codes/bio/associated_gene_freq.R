library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")
hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% rename(gene = gene_name) -> hg19
top_genes %>% group_by(gene) %>%summarize(freq = n()) -> top_genes
associated_gene_freq <- top_genes[order(top_genes$freq), ]
#View(hg19)
#View(df)
#associated_gene_freq %>% left_join(hg19, by = c("gene", "gene"))
#asscociated_gene_freq %>% select(c("gene", "freq", "gene_symbol")) -> associated_gene_freq
most_common_genes <- asscociated_gene_freq[associated_gene_freq$freq > 3, ]
View(most_common_genes)