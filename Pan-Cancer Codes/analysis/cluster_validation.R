library(cluster)
library(factoextra)
library(ggplot2)
library(ggthemes)

mat_name = "allcancer_0.001"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))
labels_df = read_tsv("data/labels.tsv")

# dist = dist(data[-1])
# dist2 = as.matrix(dist)
# dist2 = as.data.frame(dist2)
# write_tsv(dist2, file.path("matrices/sgnf_genes/", paste0(mat_name, "_dist.tsv.gz")))

clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls = clusters %>% unite(class, class1:class3, sep = ".")
cls %>% mutate(class = sub(".NA", "", class)) -> cls
cls %>% mutate(class = as.numeric(factor(class))) -> cls
rownames(cls) = cls$icgc_sample_id

sil <- silhouette(x = cls$class, distance)
summary(sil1)
fviz_silhouette(sil2, print.summary = TRUE)
  # scale_fill_brewer(palette = "Dark") +
  # scale_color_brewer(palette = "Dark") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme_minimal()
avg_sil = mean(sil[, 3])


# - -----------------------------------------------------------------------

# cls_sil = data.frame("")


# distance = read_tsv("matrices/sgnf_genes/", paste0(mat_name, "_dist.tsv.gz"))
distance = read_tsv("matrices/sgnf_genes/allcancer_0.001_dist.tsv.gz")

d = 3
clusters %>% filter(class1 == 2 & class2 == 1) -> clusters
cls = clusters %>% unite(class, class1:class3, sep = ".")
cls %>% mutate(class_2 = sub(".NA", "", class)) -> cls
clusters %>% left_join(cls) -> cls
cls %>% unite(class_2, class1:class4, sep = ".") -> cls

classes = unique(cls$class)

cls %>% 
  # arrange(class) %>% 
  mutate(class = as.numeric(factor(class))) -> c1
rownames(c1) = cls$icgc_sample_id
sil1 <- silhouette(x = c1$class, distance)
avg_sil1 = mean(sil1[, 3])

n = nrow(data)
sil_df = data.frame(class = character(), avg_sil = double())

for (c in classes) {
  if (cls %>% filter(class == c) %>% nrow() > 0.01 * n) {
  # if (cls %>% filter(class == c) %>% nrow() > 0.01 * n && grepl("NA", c) == FALSE) {
    # compute sil without applying mclust on cluster no c 
    # after applying mclust on cluster no c
    cls %>% 
      mutate(classs = ifelse(class == c, class_2, class)) %>%
      mutate(class = as.numeric(factor(classs))) -> c2
    
    cls %>% mutate(class = as.numeric(factor(class))) -> c2
    rownames(c2) = cls$icgc_sample_id
    sil2 <- silhouette(x = c2$class, distance)
    avg_sil2 = mean(sil2[, 3])
    cat("c: ", c, "  sil1: ", avg_sil1, "  sil2: ", avg_sil2, "\n")
    sil_df = rbind(sil_df, data.frame(class = c, avg_sil = avg_sil2))
  }
}

sil_df %>% mutate(class = as.character(class)) -> sil_df2

ggplot(sil_df2) + geom_bar(aes(x = class, y = avg_sil), stat = "identity", fill= "steel blue") +
  geom_abline(slope = 0, intercept = avg_sil1, color = "black", size = 1.5) + theme_calc() + scale_colour_calc() +
  ggtitle("Average silhouette")

plot_name = paste(mat_name, "(d=3) (opt)")
ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_sil.png")), width = 150, units = "mm")




