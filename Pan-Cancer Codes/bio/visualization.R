library(readr)
library(plotly)
library(factoextra)
library(RColorBrewer)

mat_name = "allcancer_0.01"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))
# clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters3.tsv")
# cls = clusters %>% unite(class, class1:class2, sep = ".") %>% 
# mutate(class = ifelse(class == "2.1", paste(class, class3, sep = "."), class)) %>% select(-class3)

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))
View(cls)
cols <- colorRampPalette(brewer.pal(9, "Set1"))
cls_plot <- function(cls) {
  ngroup = length(unique(cls$cancer_type))
  
  # plot(
  cls$class <- factor(cls$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
  ggplot(cls) + geom_histogram(aes(x = factor(class), group = cancer_type, fill = cancer_type), stat = "count") +
    scale_fill_manual(values = c("red", "palevioletred3","plum2", "turquoise1" , "mediumorchid1",
                                 "yellowgreen", "darkseagreen1", "darkgoldenrod1", "cyan4", 
                                 "blue4", "brown4","deeppink2",
                                 "chocolate1",  "darkslategrey","rosybrown",
                                   "springgreen4", "tan4",
                                 "seashell4", "black")) + xlab("Subtype") + ylab("Count") +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    guides(fill=guide_legend(title="Cancer Type")) +
      labs(title = "Figure 1.b : Cancer type distribution",
           tag = "")+
    theme(plot.tag.position = 'top',
          plot.tag = element_text(vjust = 1, hjust = 1))
    
  # )
  ggsave("bio_res/clust_histogram.pdf", width = 250, units = "mm")
  
  cls %>% group_by(class, cancer_type) %>% summarise(count = n()) -> cls_gathered
  
  plot(ggplot(cls_gathered) + geom_tile(aes(y = cancer_type, x = factor(class), fill = count)) +
         scale_fill_gradient(low = "white", high = "steelblue")) + theme_minimal()
}
cls_plot(cls)
