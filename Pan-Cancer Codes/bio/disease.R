library(readr)

path = "bio_res/gene analysis/disease"
clust_folders = list.dirs(path, recursive = FALSE, full.names = FALSE)

disease = data.frame(matrix(ncol = 4))
colnames(disease) = c('class', 'geneSet', 'description', 'pValue')
for (c in clust_folders) {
  print(c)
  dir = list.dirs(file.path(path, c), recursive = FALSE)[1]
  files = list.files(dir)
  f = files[which(startsWith(files, 'enrichment_results'))]
  if (length(f) != 0) {
    enrichment_result = read_tsv(file.path(dir, f))
    disease = rbind(disease,
                    enrichment_result %>% mutate(class = c) %>% select(class, geneSet, description, pValue)) 
  }
}
disease %>% drop_na() -> disease
disease$class = paste0("C", as.character(as.numeric(as.factor(disease$class))))

disease %>% 
  select(-description) %>% 
  mutate(bool = TRUE) %>% 
  spread(geneSet, bool, fill = FALSE) -> disease_mat

disease %>% 
  group_by(geneSet) %>% 
  summarise(total = n()) -> disease_total

# first plots -------------------------------------------------------------
disease %>% left_join(disease_total) -> disease_plot
disease_plot$class <- factor(disease_plot$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(disease_plot) + 
  geom_tile(aes(x = reorder(geneSet, total), y = class, fill = pValue), color = "grey") +
  theme_minimal() + xlab('GeneSet ID') + ylab('Subtype')+
  labs(title = "",
       subtitle = "",
       caption = "",
       tag = "Figure: Disease association with each subtype",
       fill = "pValue") +

  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + coord_fixed() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  

ggsave(file.path(paste0(path, "_tables"), "disease_geneset.pdf"), width = 600, units = "mm")

ggplot(disease_plot) + 
  geom_tile(aes(x = reorder(description, total), y = class, fill = pValue), color = "grey") +
  theme_minimal() + ggtitle('Disease heatmap') + 
  labs(title = "",
       subtitle = "",
       caption = "",
       x = "GeneSet Description",
       y = "Subtype",
       tag = "Figure: Gene set disease association with each subtype",
       fill = "p-value")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + coord_fixed()

ggsave(file.path(paste0(path, "_tables"), "disease_desc_analysis.pdf"), width = 760, height = 300,  units = "mm")



# tables ------------------------------------------------------------------

unique_diseases =
  disease %>% 
  group_by(geneSet) %>% summarise(count = n()) %>% filter(count == 1) %>% .$geneSet

write_csv(disease %>% filter(geneSet %in% unique_diseases),
          file.path(paste0(path, "_tables"), "unique_diseases.csv"))

write_csv(
  disease %>% filter(!(geneSet %in% unique_diseases)) %>% 
    spread(class, pValue, fill = 1) ,
  file.path(paste0(path, "_tables"), "common_diseases.csv"))


disease %>% left_join(disease_total) %>% spread(class, pValue, fill = 1) %>% 
  gather("class", "pValue", `1.1`:`2.2`) -> disease_plt


write_csv(disease_plot,
          file.path(paste0(path, "_tables"), "disease.csv"))

