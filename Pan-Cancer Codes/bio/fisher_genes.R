{
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pbapply)
  library(data.table)
}

# feature genes -------------------------------------------------------------------

mat_name = "allcancer_0.001"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))
labels_df = read_tsv("data/labels.tsv")

data %>% 
  gather(gene, mutation_count, ENSG00000002746:ENSG00000254709) %>% 
  filter(mutation_count != 0) -> gene_mat

#View(gene_mat)
# all genes ---------------------------------------------------------------

sample_data = read_tsv("data/sample_data.tsv")
hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) -> hg19

sample_data %>% 
  rename(gene = gene_affected) %>% 
  left_join(hg19) %>% filter(Gene_class == 'protein_coding') %>% select(-Gene_class) %>% 
  group_by(icgc_sample_id, gene) %>% 
  summarise(mutation_count = n()) -> gene_mat


# write_tsv(gene_mat, "data/temp/codings_mat.tsv.gz")

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))

total_inclass = cls %>% group_by(class) %>% summarise(inclass = n())

fisher_genes <- function(gene_mat, cls) {
  gene_mat %>% 
    left_join(cls) %>% 
    drop_na() %>% 
    select(icgc_sample_id, gene, class) %>%
    # unique() %>% 
    group_by(gene, class) %>%
    summarise(count = n()) -> genevsclass
  
  total_ingene = genevsclass %>% group_by(gene) %>% summarise(ingene = sum(count))
  
  genevsclass %>% left_join(total_inclass) %>% left_join(total_ingene) -> genevsclass
  total = total_inclass %>% .$inclass %>% sum()
  genevsclass %>% 
    mutate(not_inclass = ingene - count) %>%
    mutate(not_ingene = inclass - count) %>% 
    mutate(not_inboth = (total - inclass) - not_inclass) -> genevsclass_s
    
  # write_tsv(genevsclass_s, "data/temp/genevsclass_codings.tsv.gz")
  # genevsclass_s = read_tsv("data/temp/genevsclass_codings.tsv.gz")
  
  genevsclass_s %>% ungroup() -> genevsclass_s
  mat = matrix(nrow = 2, ncol = 2)
  geneclust_pval <- function(g, c) {
    #             class         not in class
    # gene
    # not gene
    
    row = genevsclass_s %>% filter(class == c & gene == g)
    # row = genevsclass_s[which(genevsclass_s$class == c & genevsclass_s$gene == g), ]
    mat[1, 1] = row$count
    mat[1, 2] = row$not_inclass
    mat[2, 1] = row$not_ingene
    mat[2, 2] = row$not_inboth
    t = fisher.test(mat, alternative = "greater")
    t$p.value
  }
  
  l = pbapply(genevsclass_s[, c('gene', 'class')], 1, function(x) geneclust_pval(x[1], x[2]))
  fisher_df = cbind(genevsclass_s, pval = l)
  View(fisher_df)
  
  #fisher_df %>% select(gene, class, pval) -> fisher_df
  
  fisher_df %>% spread(class, pval, fill = 1) -> fisher_spread
  
  fisher_df
}
p <- fisher_genes(gene_mat , cls)

write_tsv(fisher_genes(gene_mat, cls), "bio_res/gene analysis/fisher/total_exp_codings.tsv.gz")
#write_tsv(fisher_genes(gene_mat, cls), "bio_res/gene analysis/fisher/fisher_df_exp_codings.tsv.gz")
