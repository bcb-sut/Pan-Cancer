library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lattice)

data = read_tsv("data/sample_data.tsv.gz")

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls

hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) -> hg19

cosmic = read_csv("data/Census_alloWed.csv")

data %>% 
  select(-c(motif_3mer, count)) %>% 
  distinct() -> gene_data

write_tsv(gene_data, "data/gene_data.tsv.gz")
#gene_data = read_tsv("data/gene_data.tsv.gz")

gene_data %>% 
  left_join(cls) %>% 
  group_by(gene_affected, class) %>% 
  summarize(count = n()) -> gene_joined

cls %>% group_by(class) %>% summarise(inclass = n()) -> inclass
gene_joined %>% left_join(inclass) %>% mutate(rate = count / inclass) %>% drop_na() -> gene_joined

plt_gene_rate <- function(label, name) {
  gene_joined %>% 
    rename(gene = gene_affected) %>% 
    left_join(hg19) %>% 
    filter(Gene_class == label) -> gene_df
  # filter(gene %in% census_genes$gene_id) -> gene_df
  
  p = ggplot(gene_df) +
    geom_bar(aes(x = gene_symbol, y = rate), stat = "identity") +
    xlab("gene symbol") +
    ggtitle("Gene rate (coding)") +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.x=element_blank()) +
    facet_wrap(class ~ .)
  
  # barchart(~count|gene, group = factor(class),data=gene_df,
  #          # key = simpleKey(text = as.character(1:n_distinct(gene_df$class)),
  #                          rectangles = TRUE, points = FALSE, space = "right")
  
  # barplot(gene_df$count, main = "Gene rate in coding genes")
  # dev.copy(jpeg, filename = "bio_res/gene analysis/gene rate plots/coding.jpeg")
  # dev.off()
  
  ggsave(file.path("bio_res/gene analysis/gene rate plots", paste0(label, ".png")), 
         plot = p, width = 1200, height = 500, units = "mm")
}

cosmic %>% 
  select(Synonyms) %>% 
  rename(gene_id = Synonyms) %>% 
  separate_rows(gene_id, sep = ",") %>% 
  filter(gene_id %in% hg19$gene) -> census_genes

plt_gene_rate(c('lincRNA', 'protein_coding'), 'coding and non-coding')


