library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspectra)


######################################### pre define 2 types of data : 'data' , 'data_in_sgnf_genes'

# 'data' : a table shows number of mutation in each motif of each gene for all samples
data = read_tsv("data/sample_data.tsv")
View(data)
# feature genes : 684
sgnf_genes = read_tsv("sgnf_genes/allcancer_0.001.tsv")
View(sgnf_genes)
# a table like 'data' but only contains info about sgnf _genes in it
data_in_sgnf_gene <- data[which(data$gene_affected %in% sgnf_genes$gene_id), ]


#read clusering , the table assign each sample to a cluster
cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))
#View(cls)


hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) -> hg19

cosmic = read_csv("data/Census_alloWed.csv")

data_in_sgnf_gene %>% 
  select(-motif_3mer) %>% 
  left_join(cls) %>% 
  group_by(gene_affected, icgc_sample_id, class) %>%
  summarize(m_count = sum(count)) %>%
  ungroup() -> gene_joined

gene_joined <- gene_joined %>%
  select(-icgc_sample_id) %>%
  group_by(gene_affected, class) %>%
  summarise(count_m = sum(m_count))

View(gene_joined)

gene_joined %>% rename(gene = gene_affected) -> gene_joined
sgnf_genes %>% rename(gene = gene_id) -> sgnf_genes
merge(x=gene_joined %>% filter(class == "C1"), y = sgnf_genes, by = "gene", all.y = TRUE) -> all_gene_joined
all_gene_joined$class = "C1"
for (i in 2:17){
  merge(x=gene_joined %>% filter(class == paste0("C", i)), y = sgnf_genes, by = "gene", all.y = TRUE) -> temp
  temp$class = paste0("C", i)
  all_gene_joined <- rbind(all_gene_joined, temp)
}
all_gene_joined$count_m[is.na(all_gene_joined$count_m)] <- 0

all_gene_joined %>% 
  left_join(hg19) %>%
  filter(Gene_class == "protein_coding")-> gene_df
  #filter(gene %in% census_genes$gene_id) -> 
gene_df <- gene_df %>% filter(class %in% c("C1", "C2"))
test <- gene_df %>% filter (class == "C1") %>% left_join( gene_df %>% filter (class == "C2"), by= "gene") %>% mutate(frc = count_m.x / count_m.y)
p = ggplot(test %>% select(gene, frc))+
  geom_bar(aes(x = gene, 
               y = frc),
           stat = "identity", show.legend = FALSE)+
  xlab("gene symbol") +
  #stat_peaks(aes(x= gene_symbol, y = rate),  color = "red", geom = "text", hjust = -0.1)+
  #stat_label_peaks(span = 41, geom = "label")+
  theme(axis.text.x = element_blank())+
  labs(tag = "mutational rate",
       ylab = "fraction")+
  ylim(0, 2) +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1, size = 30), axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 30), strip.text = element_text(size = 30)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave(file.path("bio_res/gene analysis", paste0("frac_mutation_rate_C1_C2", ".pdf")), 
       plot = p, width = 600, height = 400, units = "mm")


#print(sum(test$count_m.y)/sum(test$count_m.x))
#gene_df$class <- factor(gene_df$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
gene_df$class <- factor(gene_df$class, levels = c("C1","C2"))
View(gene_df)
p = ggplot(gene_df) +
  geom_bar(aes(x = gene_symbol, 
               y = count_m, fill = class),
           stat = "identity", show.legend = FALSE)+
  xlab("gene symbol") +
  #stat_peaks(aes(x= gene_symbol, y = rate),  color = "red", geom = "text", hjust = -0.1)+
  #stat_label_peaks(span = 41, geom = "label")+
  facet_grid(class ~ .) +
  theme(axis.text.x = element_blank())+
  labs(tag = "mutational rate",
       ylab = "fraction")+
  ylim(0, 3000) +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1, size = 30), axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 30), strip.text = element_text(size = 30)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave(file.path("bio_res/gene analysis", paste0("mutation_rate_C1_C2", ".pdf")), 
       plot = p, width = 600, height = 400, units = "mm")



View(all_gene_joined)