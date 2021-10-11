library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))


motifs = read_tsv("data/sample_data.tsv")
motifs %>% 
  rename(gene = gene_affected) %>% 
  select(icgc_sample_id, motif_3mer, gene) %>% 
  left_join(cls) -> motif_joined

cls %>% 
  group_by(class) %>%
  summarise(n = n()) -> inclass
View(motif_joined)
motifrate_ingene <- function(gene_id, gene_symbol, subtype_list) {
  motif_joined %>% 
    filter(gene == gene_id) %>% 
    group_by(class, motif_3mer) %>% 
    summarise(motif_rate = n()) %>% 
    left_join(inclass) %>% 
    mutate(motif_rate_nrm = motif_rate / n) -> motif_plt
  
  # motif_joined %>% 
  #   filter(gene == gene_id) -> df
  # df %>% 
  #   group_by(class) %>%
  #   summarise(sum = n()) -> motif_sum
  # df %>% 
  #   group_by(class, motif_3mer) %>% 
  #   summarise(motif_rate = n()) %>% 
  #   left_join(motif_sum) %>% 
  #   mutate(motif_rate_nrm = motif_rate / sum) -> motif_plt
    
  ref_df = mutation_ref()
  ref_df %>% left_join(motif_plt) -> motif_plt
  allmotifs <- motif_plt$motif_3mer %>% unique()
  allmotifs <- data.frame(allmotifs) %>% rename(motif_3mer = allmotifs)
  merge(x=motif_plt %>% filter(class == subtype_list[1]), y = allmotifs, by = "motif_3mer", all.y = TRUE) -> new_motif_plt
  new_motif_plt$class = subtype_list[1]
  for (i in 2:6){
    merge(x=motif_plt %>% filter(class == subtype_list[i]), y = allmotifs, by = "motif_3mer", all.y = TRUE) -> temp
    temp$class = subtype_list[i]
    new_motif_plt <- rbind(new_motif_plt, temp)
  }
  new_motif_plt$motif_rate_nrm[is.na(new_motif_plt$motif_rate_nrm)] <- 0
  new_motif_plt$n[is.na(new_motif_plt$n)] <- 0
  new_motif_plt$motif_class[is.na(new_motif_plt$motif_class)] <- 0
  new_motif_plt$motif_rate[is.na(new_motif_plt$motif_rate)] <- 0
  
  
  
  motif_plt <- new_motif_plt
  View(motif_plt)
  motif_plt$class <- factor(motif_plt$class, levels = c("C4", "C6", "C9", "C13", "C14", "C15"))
  p = ggplot(motif_plt %>% drop_na()) + 
    geom_bar(aes(x = motif_3mer, y = motif_rate_nrm, fill = factor(motif_class)),
             stat = "identity", show.legend = FALSE) +
    #stat_peaks(span = 41, shape = 21, size = 3) +
    theme(axis.text.x = element_blank()) +
    ylab("motif rate in subtype") + xlab("3mer motif") +
    facet_wrap(class ~ ., strip.position="bottom") +
    labs(title = paste0('Figure 4a : Motif rate analysis in ', gene_symbol),
  fill = "")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1, size = 30), axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25), strip.text = element_text(size = 20))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
         axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  
  ggsave(file.path("bio_res/motif analysis/in frequent genes/", paste0("gene_", gene_symbol, ".pdf")), 
         plot = p, width = 500, height = 300,  units = "mm")
}

motifrate_ingene("ENSG00000155657", 'TTN')


# ggplot(motif_gathered %>% drop_na()) + 
#     geom_bar(aes(x = motif_3mer,
#                  y = motif_rate_nrm),
#              stat = "identity", fill = "steel blue") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     facet_grid(class ~ .)
# 
# ggsave(file.path(paste0("../plots/", mat_name), "motif_rate.png"), width = 400, height = 250, units = "mm")

# in cosmic genes  --------------------------------------------------------------------

# motifs %>% 
#   rename(gene = gene_affected) %>% 
#   select(icgc_sample_id, gene) -> gene_mat
# 
# fisher_df = fisher_genes(gene_mat)
# 
hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) %>% 
  dplyr::select(gene, gene_symbol) -> hg19
rownames(hg19) = hg19$gene

cosmic_genes = read_tsv("sgnf_genes/from_dataset/census_genes.tsv")
cosmic_genes = cosmic_genes$gene_id

for (g in cosmic_genes) {
  print(hg19[g, ]$gene_symbol)
  motifrate_ingene(g, hg19[g, ]$gene_symbol)
}


# in frequent genes -------------------------------------------------------

motifs %>% 
  rename(gene = gene_affected) %>% 
  select(icgc_sample_id, motif_3mer, gene) %>% 
  left_join(cls) -> motif_joined

motif_joined %>% 
  select(icgc_sample_id, gene) %>% 
  distinct() %>% 
  group_by(gene) %>% 
  summarize(total = n()) -> freq_genes

freq_genes %>% 
  drop_na() %>% 
  mutate(rank = order(order(total, decreasing = TRUE))) %>% 
  filter(rank <= 20) %>% .$gene -> freqs

View(freqs)
for (g in freqs) {
  print(hg19[g, ]$gene_symbol)
  motifrate_ingene(g, hg19[g, ]$gene_symbol)
}

#in common genes of subtype(n >=4)-----------------------------
# here we do the same analysis with genes that were sgnf in at leeast 4 subtypes
top_genes = read_tsv("bio_res/gene analysis/top genes/top100_exp_codings.tsv")
View(top_genes)
hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) %>% 
  dplyr::select(gene, gene_symbol) -> hg19
rownames(hg19) = hg19$gene

top_genes %>% group_by(gene) %>% summarize(freq = n()) -> temp
associated_gene_freq <- temp[order(temp$freq), ]

associated_gene_freq <- associated_gene_freq[associated_gene_freq$freq > 3, ]
associated_gene_freq_with_subtype <- associated_gene_freq %>% left_join(top_genes)
View(associated_gene_freq_with_subtype)

for (k in associated_gene_freq$gene){
  print(k)
  print(hg19[k, ]$gene_symbol)
  temp <- associated_gene_freq_with_subtype %>% filter(gene == k)
  motifrate_ingene(k, hg19[k, ]$gene_symbol, temp$clust)
}

motifrate_ingene("ENSG00000078328", 'RBFOX1')
motifrate_ingene("ENSG00000120322", 'PCDHB8')
motifrate_ingene("ENSG00000141510", 'TP53')
motifrate_ingene("ENSG00000155657", 'TTN')
motifrate_ingene("ENSG00000185507", 'IRF7')
motifrate_ingene("ENSG00000185915", 'KLHL34')
motifrate_ingene("ENSG00000205281", 'GOLGA6L10')
motifrate_ingene("ENSG00000213281", 'NRAS')
motifrate_ingene("ENSG00000185567", 'AHNAK2')
motifrate_ingene("ENSG00000181143", 'MUC16', c("C4", "C6", "C9", "C13", "C14", "C15"))



