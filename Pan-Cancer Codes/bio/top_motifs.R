library(readr)
library(dplyr)
library(pbapply)
library(tidyr)
library(ggplot2)

motifs = read_tsv("data/sample_data.tsv")

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% dplyr::select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))
View(cls)


motifs %>% 
  rename(motif = motif_3mer) %>% 
  select(icgc_sample_id, motif) %>% 
  distinct() %>% 
  left_join(cls) %>% 
  group_by(class, motif) %>% 
  summarise(count = n()) %>% 
  mutate(rank = order(order(count, decreasing = TRUE))) %>% 
  filter(rank <= 50) %>% 
  left_join(n_clusts) %>% mutate(diff = n - count) -> a

# Fisher  -----------------------------------------------------------------
sgnf_genes = read_tsv("sgnf_genes/allcancer_0.001.tsv")
top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")
View(top_genes)
View(motifs)
#View(top_genes$ gene %>% unique())
motifs %>% rename(motif = motif_3mer) %>% 
  dplyr::select(icgc_sample_id, motif, gene_affected) %>% 
  distinct() %>% 
  left_join(cls) %>% 
  drop_na() -> motifs
View(motifs)
# a table like 'data' but only contains info about sgnf _genes in it
motifs %>% inner_join(top_genes, by = c("gene_affected" = "gene", "class"= "clust")) %>% group_by(motif, class, gene_affected) %>%
  summarise(count = n()) %>% 
  filter(!grepl('N', motif)) -> motifvsclass
View(motifvsclass)
total_inmotif = motifvsclass %>% group_by(motif, gene_affected) %>% summarise(inmotif = sum(count))
total_inclass = cls %>% group_by(class) %>% summarise(inclass = n())
motifvsclass %>% left_join(total_inclass) %>% left_join(total_inmotif) -> motifvsclass
total = total_inclass %>% .$inclass %>% sum()
motifvsclass %>% 
  mutate(not_inclass = inmotif - count) %>%
  mutate(not_inmotif = inclass - count) %>% 
  mutate(not_inboth = (total - inclass) - not_inclass) -> motifvsclass_s  


motifvsclass_s %>% ungroup() -> motifvsclass_s
View(motifvsclass_s)
mat = matrix(nrow = 2, ncol = 2)
motifclust_pval <- function(m, c, g) {
  #             class         not in class
  # motif
  # not motif
  
  row = motifvsclass_s %>% filter(motif == m & class == c & gene_affected ==g)
  # row = motifvsclass_s[which(motifvsclass_s$class == c & motifvsclass_s$motif == g), ]
  mat[1, 1] = row$count
  mat[1, 2] = row$not_inclass
  mat[2, 1] = row$not_inmotif
  mat[2, 2] = row$not_inboth
  t = fisher.test(mat, alternative = "greater")
  t$p.value
}

l = pbapply(motifvsclass_s[, c('motif', 'class', 'gene_affected')], 1, function(x) motifclust_pval(x[1], x[2], x[3]))
fisher_df = cbind(motifvsclass_s, pval = l)

fisher_df %>% select(motif, class, gene_affected,pval) -> fisher_df
View(fisher_df)
fisher_df %>% 
  rename(clust = class) %>% 
  group_by(clust) %>% 
  arrange(pval) %>%
  mutate(rank = order(order(-pval, decreasing = TRUE))) %>% 
  filter(rank <= 100) -> top_motifs
View(top_motifs)
for (i in (1:17)){
  write_csv(top_motifs %>% filter(clust == paste0("C", i))  %>% select(clust,gene_affected, motif, pval), paste0("bio_res/motif analysis/", i, "top motif.csv"))
}
write_csv(top_motifs, "bio_res/motif analysis/top motif.csv")

top_motifs$gene <- paste(top_motifs$motif, top_motifs$gene_affected)


topgenes_interaction <- function(top_genes) {
  top_genes %>% 
    mutate(has_gene = TRUE) %>% 
    dplyr::select(c(has_gene, clust, gene)) %>% 
    spread(clust, has_gene, fill = FALSE) -> gene_clust
  
  cols = colnames(gene_clust)[-1]
  n = length(cols)
  interaction_df = data.frame(matrix(ncol = n, nrow = n))
  colnames(interaction_df) = cols
  rownames(interaction_df) = cols
  for (c1 in cols) {
    for (c2 in cols) {
      interaction_df[c1, c2] = sum(gene_clust[c1] & gene_clust[c2])
    }
  }
  interaction_df$c1 = rownames(interaction_df)
  interaction_df %>% gather("c2", "interactions", 1:n) %>% unique() -> plt_df
  plt_df$c1 <- factor(plt_df$c1, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
  plt_df$c2 <- factor(plt_df$c2, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
  
  
  p <- ggplot(data = plt_df, mapping = aes(x = c1, y = c2, fill = interactions))+
    geom_tile(size = 4) + geom_text(aes(label = interactions), size = 2) +
    scale_fill_gradient2(low = "lavender", high = "blueviolet",
                         name="interactions") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle =0))+
    labs(title = "Figure3b: Gene-motif interaction between each two subtypes",
         subtitle = "",
         caption = "",
         x = "Subtype", y = "Subtype",
         tag = "",
         fill = "fraction")+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))
  p
  ggsave(file.path("bio_res/motif analysis", "interactions_gene_motif_100_exp.pdf"), width = 200, units = "mm")
  
}

topgenes_interaction(top_motifs)






