{
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
}

mat_name = "allcancer_0.01"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))

clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters2.tsv")
cls = clusters %>% unite(class, class1:class2, sep = ".") %>% 
  mutate(cluster = as.numeric(as.factor(class)))

clust.centroid = function(i, data, clusters) {
  ind = (clusters == i)
  colMeans(data[ind,])
}
centers = sapply(unique(cls$cluster), clust.centroid, data[, -1], cls$cluster)

get_elbow <- function(x, cls) {
  ss <- function(x) sum(scale(x, scale = FALSE) ^ 2)
  withinss <- sapply(split(as.data.frame(x), cls$cluster), ss)
  tot.withinss <- sum(withinss) # or  resid <- x - fitted(cl); ss(resid)
  tot.withinss
}

clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters3.tsv")
depth_cluster <- function(k) { # k is actually depth
  if (k == 1) {
    cls = clusters %>% mutate(class = class1) %>% 
      mutate(cluster = as.numeric(as.factor(class)))
  }
  if (k == 2) {
    cls = clusters %>% unite(class, class1:class2, sep = ".") %>% 
      mutate(cluster = as.numeric(as.factor(class)))
  }
  if (k == 3) {
    cls = clusters %>% unite(class, class1:class2, sep = ".") %>% 
      mutate(class = ifelse(class == "1.5", paste(class, class3, sep = "."), class)) %>% 
      mutate(cluster = as.numeric(as.factor(class)))
  }
  x = data[, -1]
  get_elbow(x, cls)
}

elbow = data.frame(depth = 1:3, wss = sapply(1:3, depth_cluster))
ggplot(elbow, aes(x = depth, y = wss)) + geom_line(color = "steel blue") + 
  geom_point(color = "steel blue") + theme_bw() + ggtitle("clustering 1.5")

# - -----------------------------------------------------------------------

# clusters %>% filter(class1 == 2 & class2 == 1) -> clusters
clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters3.tsv")
cls = clusters %>% unite(class, class1:class2, sep = ".")
cls %>% mutate(class_2 = sub(".NA", "", class)) -> cls
clusters %>% left_join(cls) -> cls
cls %>% unite(class_2, class1:class3, sep = ".") -> cls

classes = unique(cls$class)

cluster_class <- function(c) {
    cls %>% 
      mutate(classs = ifelse(class == c, class_2, class)) %>%
      mutate(cluster = as.numeric(as.factor(classs))) -> c2
  x = data[, -1]
  get_elbow(x, c2)
}

elbow = data.frame(class = classes, wss = sapply(classes, cluster_class))
ggplot(elbow, aes(x = class, y = wss)) + geom_line(color = "steel blue", group = 1) + 
  geom_point(color = "steel blue") + theme_bw() + ggtitle("wss for clusters in depth 3")


