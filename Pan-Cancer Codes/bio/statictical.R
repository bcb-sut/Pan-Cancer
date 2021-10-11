library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))
View(cls)


ngroup = length(unique(cls$cancer_type))
cols <- colorRampPalette(brewer.pal(9, "Set1"))
ggplot(cls) +
  geom_histogram(aes(x = factor(class), group = cancer_type, fill = cancer_type), stat = "count") +
  scale_fill_manual(values = cols(ngroup))

cls %>% group_by(class, cancer_type) %>% summarise(count = n()) -> cls_gathered
View(cls_gathered)
cls_gathered$class <- factor(cls_gathered$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(cls_gathered, aes(y = cancer_type, x = class)) + 
  geom_tile(aes(fill = count), size=4) +
  geom_text(aes(label = count), size=3) +
  scale_fill_gradient2(low = "lavender", high = "steelblue", na.value = "0", space = "Lab") + 
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())+
  labs(title = "FigureS1: fraction of each cancer type in identified subtypes",
       subtitle = "",
       caption = "",
       x = "Subtype", y = "Cancer type",
       tag = "",
       fill = "count")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))
  

ggsave(file.path("bio_res", "cancertype_distribution_heatmap.pdf"), width = 200, units = "mm")


cls %>% group_by(class) %>% summarize(inclass = n()) -> inclass
cls_gathered %>% left_join(inclass) %>% mutate(freq = round(count / inclass, 3) * 100) -> cls_plt
cls_plt$class <- factor(cls_plt$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(cls_plt, aes(y = cancer_type, x = class)) + 
  geom_tile(aes(fill = freq), size=4) +
  geom_text(aes(label = freq), size=3) +
  scale_fill_gradient(low = "lavender", high = "blueviolet", na.value = "0") + 
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())+
  labs(title = "FigureS1: fraction of each cancer type in identified subtypes",
       subtitle = "",
       caption = "",
       x = "Subtype", y = "Cancer type",
       tag = "",
       fill = "fraction")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))

ggsave(file.path("bio_res", "cancertype_distribution_heatmap_normalized.pdf"), width = 250, units = "mm")

