library(readr)
library(dplyr)
library(tidyr)
library(pbapply)

motifs = read_tsv("data/sample_data.tsv")

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls

top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")

motifs %>% 
  rename(gene = gene_affected) %>% 
  rename(motif = motif_3mer) %>% 
  select(icgc_sample_id, motif, gene) %>% 
  distinct() %>% 
  left_join(cls) -> motif_joined

motif_joined %>% 
  filter(gene %in% top_genes$gene) %>% 
  group_by(class, gene, motif) %>% 
  summarize(count = n()) -> df  # count = number of patients mutated in that class, gene, and motif

# df %>% mutate(hasRow = TRUE) %>% spread(class, hasRow, fill = FALSE) %>% View()

df %>% group_by(gene, motif) %>% summarise(n_class = n()) %>% filter(n_class > 1) -> common_gm
classes = unique(cls$class)

motif_joined %>% 
  filter(gene %in% top_genes$gene) %>% 
  select(class, gene, icgc_sample_id) %>% 
  distinct() %>% 
  group_by(class, gene) %>%
  summarise(in_gene = n()) -> gene_df

common_gm %>% left_join(df) %>% rename(c1 = class) %>% rename(count1 = count) %>% 
  left_join(df) %>% rename(c2 = class) %>% rename(count2 = count) %>% 
  filter(c1 < c2) %>% 
  left_join(gene_df, by = c("c1" = "class", "gene")) %>% rename(in_gene1 = in_gene) %>% 
  left_join(gene_df, by = c("c2" = "class", "gene")) %>% rename(in_gene2 = in_gene) %>% 
  mutate(ncount1 = in_gene1 - count1) %>% 
  mutate(ncount2 = in_gene2 - count2) %>%
  select(c1, c2, gene, motif, in_gene1, count1, ncount1, in_gene2, count2, ncount2) %>% 
  distinct() -> gm_df

gm_df = gm_df %>% ungroup()
mat = matrix(nrow = 2, ncol = 2)
motifass_pval <- function(count1, ncount1, count2, ncount2) {
  
  #     in motif  not in motif
  #     count1    ncount1
  #     count2    ncount2
  
  # row = gm_df %>% filter(gene == g & motif == m & c1 == clust1 & c2 == clust2)
  mat[1, 1] = count1
  mat[1, 2] = ncount1
  mat[2, 1] = count2
  mat[2, 2] = ncount2
  t = fisher.test(mat, alternative = "two.sided")
  t$p.value
}

l = pbapply(gm_df[, c("count1", "ncount1", "count2", "ncount2")], 1, 
            function(x) motifass_pval(x[1], x[2], x[3], x[4]))
fisher_df = cbind(gm_df, pval = l)


write_csv(fisher_df, "bio_res/motif analysis/pairwise/motif_fisher.csv")



