{
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
}

labels_df = read_tsv("data/labels.tsv")

labels_df %>% 
  group_by(cancer_type) %>% 
  summarise(n = n()) -> cancer_count

clusters = clust %>% left_join(labels_df) %>% mutate(class1 = class) %>% select(-class)
# write_tsv(clusters, "method_res/model_based/allcancer_0.001/clusters1.tsv")
# 
# clusters %>% filter(class1 == 1) %>%  .$icgc_sample_id -> c1
# clust1 = cbind(icgc_sample_id = c1, class = mc1$classification) %>% as.data.frame()
# clusters %>% filter(class1 == 2) %>%  .$icgc_sample_id -> c2
# clust2 = cbind(icgc_sample_id = c2, class = mc2$classification) %>% as.data.frame()
# clust1 = rbind(clust1, clust2)
# 
# clusters = clusters %>% left_join(clust1) %>% drop_na() %>% rename(class2 = class)
# write_tsv(clusters, "method_res/model_based/allcancer_0.001/clusters2.tsv")
# 

cls = clusters %>% unite(class, class1:class2, sep = ".")
cls_plot(cls, "allcancer_0.001 (d=2)")

cls_plot <- function(cls, plot_name) {
  
  ngroup = length(unique(cls$cancer_type))
  p1 <- ggplot(cls) + geom_histogram(aes(x = factor(class), fill = cancer_type), stat = "count") + 
    scale_fill_manual(values = cols(ngroup)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("class") + ggtitle(plot_name)
  plot(p1)
  ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_bar.png")), width = 250, units = "mm")
  cls %>% group_by(class, cancer_type) %>% summarise(count = n()) %>% left_join(cancer_count) -> cls_gathered
  
  p2 <- ggplot(cls_gathered) + geom_tile(aes(y = cancer_type, x = factor(class), fill = count)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient(low = "white", high = "steelblue") + ggtitle(plot_name) + xlab("class")
  plot(p2)
  ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_heatmap.png")), width = 200, units = "mm")
  
  p3 <- ggplot(cls_gathered) + geom_tile(aes(y = cancer_type, x = factor(class), fill = (count / n))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient(low = "white", high = "steelblue") + ggtitle(plot_name) + xlab("class")
  plot(p3)
  ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_heatmap_freq.png")), width = 200, units = "mm")
}

subclust <- function(labels) {
  data %>% filter(icgc_sample_id %in% labels) -> data1
  rownames(data1) = data1$icgc_sample_id
  mc1 = Mclust(data1[, -1])
  clust = data.frame(cbind(icgc_sample_id = data1$icgc_sample_id, class = mc1$classification), stringsAsFactors = FALSE)
  clust
}

clusters %>% filter(class1 == 1) %>% .$icgc_sample_id -> c1
clust2 = subclust(c1)
clusters %>% filter(class1 == 2) %>% .$icgc_sample_id -> c2
clust22 = subclust(c2)
clust1 = rbind(clust2, clust22)
clusters = clusters %>% left_join(clust1) %>% drop_na() %>% rename(class2 = class)

# write_tsv(clusters, "method_res/model_based/cancus/clusters3.tsv")
clusters = read_tsv("method_res/model_based/cancus/clusters2.tsv")

mat_name = "cosmic"
cls = clusters %>% unite(class, class1:class3, sep = ".")
cls %>% mutate(class = sub(".NA", "", class)) -> cls
cls_plot(cls, "cosmic (d=3)")

clusters %>% filter(class1 == 2 & class2 == 1) %>% .$icgc_sample_id -> c1
clust2 = subclust(c1)
clusters %>% filter(class1 == 2 & class2 == 2) %>% .$icgc_sample_id -> c2
clust22 = subclust(c2)
clust1 = rbind(clust2, clust22)
clusters = clusters %>% left_join(clust1) %>% rename(class3 = class)
clusters = read_tsv("method_res/model_based/cancus/clusters3.tsv")

