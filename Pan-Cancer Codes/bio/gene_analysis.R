library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

fisher_df = read_tsv("bio_res/gene analysis/fisher/total_exp_codings_complete.tsv")
#View(fisher_df)
top_n <- function(n) {
  fisher_df %>% 
    rename(clust = class) %>% 
    group_by(clust) %>% 
    arrange(pval) %>%
    mutate(rank = order(order(-pval, decreasing = TRUE))) %>% 
    filter(rank <= n) -> top_genes
  
  top_genes
}

top_genes = top_n(100)
View(top_genes)
write_tsv(top_genes, "bio_res/gene analysis/top genes/top100_exp_codings_complete.tsv")


# Top genes for each cluster ----------------------------------------------

top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")
clusters = unique(top_genes$clust)
for (c in clusters) {
  d = dir.create(paste("bio_res/gene analysis/top genes/exp", c, sep = "/"))
  write_tsv(top_genes %>% filter(clust == c),
            # %>% ungroup() %>% select(gene),
            file.path(paste("bio_res/gene analysis/top genes/exp", c, sep = "/"), "top_genes.tsv"))
}

# Interactions heatmap ----------------------------------------------------

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
    geom_tile(size = 5) + geom_text(aes(label = interactions), size = 3) +
    scale_fill_gradient2(low = "lavender", high = "blueviolet",
                         name="interactions") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 0))+
    labs(title = "Figure3.a: Gene interaction between each two subtypes",
         subtitle = "",
         caption = "",
         x = "Subtype", y = "Subtype",
         tag = "",
         fill = "fraction")+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))
  p
}

top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")

topgenes_interaction(top_genes)

ggsave(file.path("bio_res/gene analysis/top genes/", "interactions_100_exp.pdf"), width = 200, units = "mm")


# ranking -----------------------------------------------------------------

top_genes %>% 
  select(gene, clust) %>% 
  group_by(gene) %>% 
  summarize(appeared = n()) %>% 
  arrange(-appeared) -> common_genes

  # left_join(top_genes) -> common_genes

hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) %>% 
  dplyr::select(gene, gene_symbol) -> hg19


View(hg19)
View(common_genes)
common_genes %>% 
  left_join(hg19) -> common_genes



for (i in 2:16) {
  motifrate_ingene(
  common_genes[i, ]$gene,
  common_genes[i, ]$gene_symbol)
}

#################### tables for each subtype containing pvalue and name of gene ###########
top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")
top_genes <- top_genes %>% left_join(hg19)
View(top_genes)
for(i in 1:17){
  subtype_i  <- top_genes %>% filter(clust == paste0("C", i)) %>% select(c(gene_symbol, pval))
  write_csv(subtype_i, path = paste0("bio_res/gene analysis/top genes/C", i, ".csv"))
  
  
}



########################### fraction of top00 in sgnf each subtype separately #################
sgnf_genes = read_tsv("sgnf_genes/allcancer_0.001.tsv")
sgnf_genes %>% rename(gene = gene_id) -> sgnf_genes
top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")
View(sgnf_genes)
top_genes %>% group_by(clust) %>% inner_join(sgnf_genes) %>% summarise(count = n()) -> fract
View(fract)
