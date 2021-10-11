library(readr)
library(ggplot2)

path = "bio_res/gene analysis/pathway"
clust_folders = list.dirs(path, recursive = FALSE, full.names = FALSE)

pathway = data.frame(matrix(ncol = 4))
colnames(pathway) = c('class', 'geneSet', 'description', 'pValue')
for (c in clust_folders) {
  print(c)
  dir = list.dirs(file.path(path, c), recursive = FALSE)[1]
  files = list.files(dir)
  f = files[which(startsWith(files, 'enrichment_results'))]
  if (length(f) != 0) {
    enrichment_result = read_tsv(file.path(dir, f))
    pathway = rbind(pathway,
                    enrichment_result %>% mutate(class = c) %>% select(class, geneSet, description, pValue)) 
  }
}
pathway %>% drop_na() -> pathway

pathway %>% 
  select(-description) %>% 
  mutate(bool = TRUE) %>% 
  spread(geneSet, bool, fill = FALSE) -> pathway_mat

pathway %>% 
  group_by(geneSet) %>% 
  summarise(total = n()) -> pathway_total


# first plots -------------------------------------------------------------
pathway %>% left_join(pathway_total) -> pathway_plot
pathway_plot$class <- factor(pathway_plot$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(pathway_plot) + 
  geom_tile(aes(x = reorder(geneSet, total), y = class, fill = pValue), color = "grey") +
  theme_minimal() +
  labs(title = "",
       subtitle = "",
       caption = "",
       x = "GeneSet ID",
       y = "Subtype",
       tag = "Figure: Gene set Pathway association with each subtype",
       fill = "pValue") +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + coord_fixed()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggsave(file.path(paste0(path, "_tables"), "pathway_geneset.pdf"), width = 200, units = "mm")

ggplot(pathway_plot) + 
  geom_tile(aes(x = reorder(description, total), y = class, fill = pValue), color = "grey") +
  theme_minimal() + 
  labs(title = "",
       subtitle = "",
       caption = "",
       x = "GeneSet Description",
       y = "Subtype",
       tag = "Figure: Gene set Pathway association with each subtype",
       fill = "pValue") +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + coord_fixed()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggsave(file.path(paste0(path, "_tables"), "pathway_desc.pdf"), width = 200, units = "mm")



# tables ------------------------------------------------------------------
View(pathway_plot)
unique_pathways =
  pathway %>% 
  group_by(geneSet) %>% summarise(count = n()) %>% filter(count == 1) %>% .$geneSet

write_csv(pathway ,
          file.path(paste0(path, "_tables"), "pathways.csv"))

write_csv(
  pathway %>% filter(!(geneSet %in% unique_pathways)) %>% 
  spread(class, pValue, fill = 1) ,
  file.path(paste0(path, "_tables"), "common_pathways.csv"))


#pathway %>% left_join(pathway_total) %>% spread(class, pValue, fill = 1) %>% 
#  gather("class", "pValue", `1.1`:`2.2`) -> pathway_plt

















path = "bio_res/gene analysis/pathway"


temp <- read_csv("bio_res/gene analysis/pathway/KEGG_enrichr.csv")
temp$class = factor(temp$class , levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(temp) + 
  geom_tile(aes(x = Term, y = class, fill = P.value), color = "grey") +
  theme_minimal() +
  labs(title = "Figure 4.a : pathway analysis",
       subtitle = "",
       caption = "",
       x = "GeneSet Description",
       y = "Subtype",
       tag = "",
       fill = "pValue") +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank()) + coord_fixed() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))


ggsave(file.path(paste0(path, "_tables"), "KEGG_pathway.pdf"), width = 400, units = "mm")

















