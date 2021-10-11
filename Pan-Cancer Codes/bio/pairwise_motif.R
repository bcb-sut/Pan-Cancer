library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))


# motif dfs ---------------------------------------------------------------

motifs = read_tsv("data/sample_data.tsv")

motifs %>% 
  rename(gene = gene_affected) %>% 
  select(icgc_sample_id, motif_3mer, gene) %>% 
  left_join(cls) -> motif_joined

# motif_joined %>% 
# group_by(class, motif_3mer) %>% 
# summarise(motif_rate = n()) -> motif_gathered

motif_ingene_plt <- function(genes, c1, c2, path) {
  print(genes)
  print(c1)
  print(c2)
  n_genes = length(genes)
  print(n_genes)
  motif_joined %>% 
    filter(class == c1 | class == c2) %>% 
    filter(gene %in% genes) -> df
  if(dim(df) == 0){
    return("null df")
  }
  df %>% 
    group_by(class) %>%
    summarise(sum = n()) -> motif_sum
  
  df %>% 
    group_by(class, motif_3mer) %>% 
    summarise(motif_rate = n()) %>% 
    left_join(motif_sum) %>% 
    mutate(motif_rate_nrm = motif_rate / sum) -> motif_plt
  
  ref_df = mutation_ref()
  ref_df %>% left_join(motif_plt) %>% drop_na() -> motif_plt
  View(motif_plt)
  motif_plt$class <- factor(motif_plt$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
  p_bar = ggplot(motif_plt) + 
    geom_bar(aes(x = motif_3mer, y = motif_rate_nrm, fill = factor(motif_class)),
             stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 5)) +
    ylab("fraction") + 
    xlab("3mer motif")+
    labs(title = paste0("Figure : Motif rate among common genes of ", c1, " & ", c2, "(n = ", length(genes), ")")) +
    theme(plot.tag.position = 'top',
          plot.tag = element_text(vjust = 1, hjust = 1, size = 15), axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
    facet_grid(class ~ .)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  # p_bar2 = ggplot(motif_plt) + 
  #   geom_bar(aes(x = motif_3mer, y = motif_rate_nrm, fill = class),
  #            stat = "identity", position = "dodge") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   ylab("motif rate in cluster") + 
  #   xlab("3mer motif")
  
  p_line = ggplot(motif_plt) + 
    geom_line(aes(x = motif_3mer, y = motif_rate_nrm, color = class, group = class)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("fraction") + 
    xlab("3mer motif")+
    labs(title = paste0("Figure : Motif rate among common genes of ", c1, " & ", c2)) +
    theme(plot.tag.position = 'top')
  
  
  ggsave(file.path(path, paste0(c1, " & ", c2, "common_genes_bar.pdf")), plot = p_bar, width = 200, units = "mm")
  # ggsave(file.path(path, "common_genes_bar2.png"), plot = p_bar2, width = 200, units = "mm")
  #ggsave(file.path(path, "common_genes_line.png"), plot = p_line, width = 200, units = "mm")
}

motif_singlegene_plt <- function(g, symbol, c1, c2, path) {

  motif_joined %>% 
    filter(class == c1 | class == c2) %>% 
    filter(gene == g) -> df
  
  df %>% 
    group_by(class) %>%
    summarise(sum = n()) -> motif_sum
  
  df %>% 
    group_by(class, motif_3mer) %>% 
    summarise(motif_rate = n()) %>% 
    left_join(motif_sum) %>% 
    mutate(motif_rate_nrm = motif_rate / sum) -> motif_plt
  
  ref_df = mutation_ref()
  ref_df %>% left_join(motif_plt) %>% drop_na() -> motif_plt
  if(symbol =="PCDHB3"){
  View(motif_plt)
  }
  p_bar = ggplot(motif_plt) + 
    geom_bar(aes(x = motif_3mer, y = motif_rate_nrm, fill = factor(motif_class)),
             stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("motif rate in cluster") + 
    xlab("3mer motif") +
    ggtitle(paste0('motif rate in gene ', symbol)) +
    facet_grid(class ~ .)
  
  # p_bar2 = ggplot(motif_plt) + 
  #   geom_bar(aes(x = motif_3mer, y = motif_rate_nrm, fill = class),
  #            stat = "identity", position = "dodge") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   ylab("motif rate in cluster") + 
  #   xlab("3mer motif")
  
  p_line = ggplot(motif_plt) + 
    geom_line(aes(x = motif_3mer, y = motif_rate_nrm, color = class, group = class)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("motif rate in cluster") + 
    ggtitle(paste0('motif rate in gene ', symbol)) +
    xlab("3mer motif")
  
  ggsave(file.path(path, paste0(symbol, c1, " & ",c2, "_bar.pdf")), plot = p_bar, width = 200, units = "mm")
  # ggsave(file.path(path, "common_genes_bar2.png"), plot = p_bar2, width = 200, units = "mm")
  #ggsave(file.path(path, paste0(symbol, "_line.png")), plot = p_line, width = 200, units = "mm")
}

#   -----------------------------------------------------------------------

hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) %>% 
  dplyr::select(gene, gene_symbol) -> hg19

rownames(hg19) = hg19$gene

plot_interactions <- function(c1, c2, genes) {
  
  path = "bio_res/motif analysis/pairwise/"
#  path = file.path("bio_res/motif analysis/pairwise_singlegene/", paste(c1, c2, sep = '_'))
#  if (!dir.exists(path)) {
#    dir.create(path)
#  }
  #for (g in genes) {
  #  gene_symbol = hg19[g, ]$gene_symbol
  #  motif_singlegene_plt(g, gene_symbol, c1, c2, path)
  #}
  
  # path = file.path("bio_res/motif analysis/pairwise/", paste(c1, c2, sep = '_'))
  # if (!dir.exists(path)) {
  #   dir.create(path)
  # }
   motif_ingene_plt(genes, c1, c2, path)
}


#   -----------------------------------------------------------------------

top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")

top_genes %>% 
  mutate(has_gene = TRUE) %>% 
  select(gene, has_gene, clust) %>% 
  spread(clust, has_gene, fill = FALSE) -> gene_clust
View(gene_clust)


cols = colnames(gene_clust)[-1]
n = length(cols)
for (i in 1:n) {
  for (j in (i+1):n) {
    c1 = cols[i]
    c2 = cols[j]
    if (length(gene_clust$gene[(gene_clust[c1] & gene_clust[c2])]) > 0) {
      plot_interactions(c1, c2, gene_clust$gene[(gene_clust[c1] & gene_clust[c2])])
    }
  }
}
c1 = "C3"
c2 = "C7"

plot_interactions(c1, c2, gene_clust$gene[(gene_clust[c1] & gene_clust[c2])])


