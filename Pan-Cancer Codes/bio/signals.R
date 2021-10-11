{
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
}

mat_name = "allcancer_0.001"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))
labels_df = read_tsv("data/labels.tsv")

data %>% 
  gather(gene, mutation_count, ENSG00000002746:ENSG00000254709) %>% 
  filter(mutation_count != 0) -> gene_mat

# clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters3.tsv")
# cls = clusters %>% unite(class, class1:class2, sep = ".") %>% 
#   mutate(class = ifelse(class == "2.1", paste(class, class3, sep = "."), class)) %>% select(-class3)

# clusters %>% mutate(class3 = ifelse(class1 == 1 & class2 == 5, class3, NA)) -> clusters
# cls = clusters %>% unite(class, class1:class2, sep = ".")
# cls %>% mutate(class = sub(".NA", "", class)) -> cls

gene_mat %>% 
  left_join(cls) %>% 
  drop_na() -> gene_mat
# write_tsv(gene_mat, "method_res/model_based/allcancer_0.001/genemat2.tsv.gz")

# gene_mat %>% 
  # group_by(class, gene) %>% 
  # summarise(mutation_count = sum(mutation_count)) -> gene_mat # TODO: CHANGE

gene_mat %>% 
  group_by(class, gene) %>% 
  summarise(mutation_count = n()) -> gene_mat


gene_mat %>% group_by(gene) %>% summarise(sumInGene = sum(mutation_count)) -> total_gene
gene_mat %>% group_by(class) %>% summarise(sumInClass = sum(mutation_count)) -> total_class

gene_mat %>% 
  group_by(class, gene) %>% 
  summarise(mutation_count = sum(mutation_count)) %>% 
  ungroup() %>% 
  group_by(class) %>% 
  arrange(mutation_count) -> top_genes
  # mutate(rank = order(order(mutation_count, decreasing = TRUE))) %>% 
  # filter(rank <= 50) -> top_genes
  # top_n(50) -> top_genes


# write_tsv(top_genes, "method_res/model_based/allcancer_0.001/top_sys.tsv")

hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) %>% 
  dplyr::select(gene, gene_symbol) -> hg19

top_genes %>% select(gene, class) %>% group_by(gene) %>% summarise(numberOfClasses = n()) -> no_classes

# top_genes %>% left_join(hg19) %>% dplyr::select(-mutation_count) -> top_withlabels
# top_withlabels %>% group_by(class) %>% arrange(rank) -> top_withlabels

# write_tsv(top_withlabels, "method_res/model_based/allcancer_0.001/top3_2_withlabels.tsv")

# cosmic_genes = read_tsv("sgnf_genes/from_dataset/census_genes.tsv")
# cosmic_genes = cosmic_genes$gene_id
# 
# top_genes %>% mutate(cosmic = ifelse(gene %in% cosmic_genes, TRUE, FALSE)) -> top_genes

# top_genes %>% filter(gene %in% cosmic_genes) %>% mutate(cosmic = TRUE) -> top_genes_cosmic
# top_genes %>% filter(!(gene %in% cosmic_genes)) %>% mutate(cosmic = FALSE) -> top_genes_cosmicless
# top_genes = rbind(top_genes_cosmic, 
#                   # data.frame(class = "1.1",gene = "------", mutation_count = 0, rank = 0),
#                   top_genes_cosmicless)
top_genes %>% ungroup() %>% 
  left_join(no_classes) %>% left_join(hg19) %>% 
  mutate(class = as.factor(class)) -> top_genes

plot_name = "top sgnf genes (d=2)"


# tile --------------------------------------------------------------------

p1 <- ggplot(top_genes) + geom_tile(aes(y = reorder(gene_symbol, numberOfClasses), x = factor(class), fill = mutation_count)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient(low = "steelblue", high = "red") + 
  ggtitle(plot_name) + xlab("class") + ylab("") + labs(fill = "#mutations")
plot(p1)
ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, ".png")), height = 500, units = "mm")

p2 <- ggplot(top_genes %>% left_join(total_gene)) + 
  geom_tile(aes(y = reorder(gene_symbol, numberOfClasses), x = factor(class), fill = (mutation_count / sumInGene))) +
  scale_fill_gradient(low = "steelblue", high = "red") +
  ggtitle(paste(plot_name, "(normalized in gene)")) + xlab("class") + ylab("") + labs(fill = "#mutation")
plot(p2)
ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_inGene.png")), height = 500, units = "mm")

p3 <- ggplot(top_genes %>% left_join(total_class)) + 
  geom_tile(aes(y = reorder(gene_symbol, numberOfClasses), x = factor(class), fill = (mutation_count / sumInClass))) +
  scale_fill_gradient(low = "steelblue", high = "red") +
  ggtitle(paste(plot_name, "(normalized in class)")) + xlab("class") + ylab("") + labs(fill = "#mutation")
plot(p3)
ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_inClass.png")), height = 500, units = "mm")

p4 <- ggplot(top_genes %>% left_join(total_class) %>% left_join(total_gene)) + 
  geom_tile(aes(y = reorder(gene_symbol, numberOfClasses), x = factor(class), 
                fill = ((mutation_count / sumInClass) / sumInGene))) +
  scale_fill_gradient2(low = "steelblue", mid = "red", high = "white") +
  ggtitle(paste(plot_name, "(normalized in class and gene)")) + xlab("class") + ylab("") + labs(fill = "#mutation")
plot(p4)
ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_nrm.png")), height = 500, units = "mm")


classes = unique(top_genes$class)

for (c in classes) {
  top_genes %>% filter(class == c) -> top1
  top1 %>% left_join(hg19) %>% View()
}



# bar plots ---------------------------------------------------------------

ggplot(top_genes %>% left_join(total_class)) + 
  geom_bar(aes(x = reorder(reorder(gene_symbol, -numberOfClasses), -cosmic), 
               y = mutation_count, fill = cosmic), 
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(class ~ .)
  
  # scale_fill_gradient(low = "steelblue", high = "red") + 
  # ggtitle(plot_name) + xlab("class") + ylab("") + labs(fill = "#mutations")

plot_name = paste(mat_name, "(d=2) (opt)")
ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_gene_signal_2.png")), width = 700, units = "mm")


# -------
top_genes %>% left_join(total_class) %>% left_join(total_gene) -> top_genes_plt

ggplot(top_genes_plt) +
  geom_line(aes(x = reorder(gene_symbol, -numberOfClasses), 
                y = mutation_count / sumInClass, group = class, color = class)) +
  xlab("gene symbol") +
  ggtitle("Gene rate") +
  # ggtitle("Gene rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave(file.path(paste0("../plots/", mat_name), "gene_rate_exp_nrmclassgene.png"), width = 800, units = "mm")
ggsave(file.path(paste0("../plots/", "signals/gene_rate"), "gene_rate_exp_frq.png"), width = 1200, units = "mm")

#   -----------------------------------------------------------------------

motifs = read_tsv("data/motifs.tsv.gz")
cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls

motifs %>% 
  # filter(gene_affected == "ENSG00000133703") %>%
  select(icgc_sample_id, motif_3mer) %>% 
  left_join(cls) -> motif_joined

motif_joined %>% 
  group_by(class, motif_3mer) %>% 
  summarise(motif_rate = n()) -> motif_gathered

motif_gathered %>% 
  group_by(class) %>% 
  summarise(sum = sum(motif_rate)) -> motif_sum

motif_gathered %>% 
  left_join(motif_sum) %>% 
  mutate(motif_rate_nrm = motif_rate / sum) %>% 
  drop_na() -> motif_gathered

ggplot(motif_gathered %>% drop_na()) +
  geom_line(aes(x = motif_3mer, y = motif_rate_nrm, group = class, color = class)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggplot(motif_gathered) + 
#   geom_bar(aes(x = motif_3mer, 
#                y = motif_rate), 
#            stat = "identity", fill = "steel blue") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   facet_grid(class ~ .)

ggsave(file.path(paste0("../plots/", mat_name), "gene_ENSG00000133703.png"), width = 400, units = "mm")

