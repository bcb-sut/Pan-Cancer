library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls

hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) -> hg19

cosmic = read_csv("data/Census_alloWed.csv")

gene_data = read_tsv("data/gene_data.tsv.gz")


cancers = unique(cls$class)

gene_joined %>% 
  rename(gene = gene_affected) %>% 
  left_join(hg19) %>% 
  filter(Gene_class == 'lincRNA') -> lnc_df
gene_joined %>% 
  rename(gene = gene_affected) %>% 
  left_join(hg19) %>% 
  filter(Gene_class == 'protein_coding') -> coding_df
View(coding_df)

dif_func <- function(label, c1, c2) {
  #if (label == 'lincRNA') {
  #  gene_df = lnc_df
  #}
  #else {
    gene_df = gene_joined
    gene_df %>% rename(gene = gene_affected) -> gene_df
  #}
  gene_df %>% select(gene, class, rate) -> gene_df
  #gene_df %>% group_by(class) %>% summarise(sum = sum(count)) -> gene_sum
  #gene_df %>% left_join(gene_sum) %>% mutate(rate = count / sum) -> gene_df
  
  
  gene_df %>% filter(class == c1) %>% rename(C1 = rate) %>% select(-class) %>% 
    left_join(gene_df %>% filter(class == c2) %>% rename(C2 = rate) %>% select(-class)) %>% 
    mutate(dif = C2 - C1) -> diff_df
  diff_df %>% gather("col", "val", 2:4) %>% left_join(hg19) -> plt_df
  plt_df$col[plt_df$col == "C1"] <- c1
  plt_df$col[plt_df$col == "C2"] <- c2
  View(plt_df)
  p = ggplot(plt_df) +
    geom_bar(aes(x = gene_symbol, y = val), stat = "identity") +
    xlab("gene symbol") +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.x=element_blank()) +
    facet_grid(col ~ .)+
    labs(tag = paste0("Figure : Coding gene rate  + Diffrential analysis for ", c1, " & ", c2, " subtypes"),
  fill = "count",
  y = "rate")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))
    
  
  ggsave(file.path("bio_res/gene analysis/differential analysis", paste0(label, " ", c1, " & ", c2, ".pdf")), 
         plot = p, 
         # width = 200, height = 500, units = "mm"
  )
}

n = length(cancers)
for (i in 1:n) {
  for (j in (i+1):n) {
    c1 = cancers[i]
    c2 = cancers[j]
    dif_func('protein_coding', c1, c2)
    #dif_func('lincRNA', c1, c2)
  }
}


# old ---------------------------------------------------------------------

ggplot(df2) + geom_line(aes(x = gene, y = val)) +
  facet_wrap(col ~ .)

df = diff_df
par(mfrow = c(3,1))
barplot(df$count1)
barplot(df$count2)
barplot(df$diff)
