library(readr)

path = "bio_res/gene analysis/geneontology"
clust_folders = list.dirs(path, recursive = FALSE, full.names = FALSE)

geneontology = data.frame(matrix(ncol = 4))
colnames(geneontology) = c('class', 'geneSet', 'description', 'pValue')
for (c in clust_folders) {
  dir = list.dirs(file.path(path, c), recursive = FALSE)[1]
  files = list.files(dir)
  f = files[which(startsWith(files, 'enrichment_results'))]
  if (length(f) != 0) {
  enrichment_result = read_tsv(file.path(dir, f))
  geneontology = rbind(geneontology,
                       enrichment_result %>% mutate(class = c) %>% select(class, geneSet, description, pValue)) 
  }
}
geneontology %>% drop_na() -> geneontology

geneontology %>% 
  select(-description) %>% 
  mutate(bool = TRUE) %>% 
  spread(geneSet, bool, fill = FALSE) -> geneontology_mat

geneontology %>% 
  group_by(geneSet) %>% 
  summarise(total = n()) -> geneontology_total
geneontology_plot <- geneontology %>% left_join(geneontology_total)
geneontology_plot$class <- factor(geneontology_plot$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(geneontology_plot)+ 
  geom_tile(aes(x = reorder(geneSet, total), y = class, fill = pValue), color = "grey") +
  theme_minimal() + 
  labs(title = "",
       subtitle = "",
       caption = "",
       x = "GeneSet ID",
       y = "Subtype",
       tag = "Figure: Gene set Ontology association with each subtype",
       fill = "pValue") +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(r=10, l = 0, t = 0, b = 0)), panel.grid = element_blank()) + coord_fixed() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggsave(file.path(paste0(path, "_tables"), "geneontology_geneset.pdf"), width = 600, units = "mm")

ggplot(geneontology_plot) + 
  geom_tile(aes(x = reorder(description, total), y = class, fill = pValue), color = "grey") +
  theme_minimal() +
  labs(title = "",
       subtitle = "",
       caption = "",
       x = "GeneSet Description",
       y = "Subtype",
       tag = "Figure: Gene set Ontology association with each subtype",
       fill = "pValue") +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank()) + coord_fixed() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggsave(file.path(paste0(path, "_tables"), "geneontology_desc.pdf"), width = 400, units = "mm")



attach(mtcars)
par(mfrow=c(1,3))
tempBP <- read_csv("bio_res/gene analysis/geneontology/BP_enrichr.csv")
tempMF <- read_csv("bio_res/gene analysis/geneontology/MF_enrichr.csv")
tempCC <- read_csv("bio_res/gene analysis/geneontology/CC_enrichr.csv")
temp <- read_csv("bio_res/gene analysis/geneontology/enrichr.csv")

temp %>% mutate(type = 0) -> temp
temp$type[temp$Term %in% tempBP$Term] <- "Biological_Process"
temp$type[temp$Term %in% tempMF$Term] <- "Molecular_Function"
temp$type[temp$Term %in% tempCC$Term] <- "Cellular_Component"


temp %>% mutate(color = 0) -> temp
temp$color[temp$Term %in% tempBP$Term] <- "#FF0000"
temp$color[temp$Term %in% tempMF$Term] <- "#1A80C4"
temp$color[temp$Term %in% tempCC$Term] <- "#008000"

#View(temp)

#sapply(strsplit(temp$Term, "("), "[", 1)
temp$Term <- read.table(text = temp$Term, sep = "(", as.is = TRUE)$V1

temp$class = factor(temp$class , levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
colors <- temp$color[order(temp$Term)]
write_csv(temp %>% select(Term, Adjusted.P.value, class, type), paste0(path, "_tables/geneontologys.csv"))
ggplot(temp) + 
  geom_tile(aes(x = Term, y = class, fill = P.value), color = "grey") +
  theme_minimal() +
  labs(title = "Figure 4)Gene ontology analysis",
       subtitle = "",
       caption = "",
       x = "GeneSet Description",
       y = "Subtype",
       tag = "",
       fill = "pValue") +
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = colors), panel.grid = element_blank()) + coord_fixed() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))


ggsave(file.path(paste0(path, "_tables"), "BP_geneontology.pdf"), width = 400, units = "mm")





# tables ------------------------------------------------------------------

unique_geneontologys =
  geneontology %>% 
  group_by(geneSet) %>% summarise(count = n()) %>% filter(count == 1) %>% .$geneSet

write_csv(geneontology_plot,
          file.path(paste0(path, "_tables"), "geneontologys.csv"))

write_csv(
  geneontology %>% filter(!(geneSet %in% unique_geneontologys)) %>% 
    spread(class, pValue, fill = 1) ,
  file.path(paste0(path, "_tables"), "common_geneontologys.csv"))


geneontology %>% left_join(geneontology_total) %>% spread(class, pValue, fill = 1) %>% 
  gather("class", "pValue", `1.1`:`2.2`) -> geneontology_plt


