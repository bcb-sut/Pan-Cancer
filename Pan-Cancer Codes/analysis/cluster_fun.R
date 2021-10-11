{
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
}


# - -----------------------------------------------------------------------


mat_name = "allcancer_0.001"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))
labels_df = read_tsv("data/labels.tsv")

subclust <- function(labels) {
  data %>% filter(icgc_sample_id %in% labels) -> data1
  rownames(data1) = data1$icgc_sample_id
  mc1 = Mclust(data1[, -1])
  clust = data.frame(cbind(icgc_sample_id = data1$icgc_sample_id, class = mc1$classification), stringsAsFactors = FALSE)
  clust
}

clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")

# hclust(clusters)

cls = clusters %>% unite(class, class1:class4, sep = ".")
cls %>% mutate(class = sub(".NA", "", class)) -> cls

classes = unique(cls$class)

d = 5
d_clust = data.frame(row.names = c("icgc_sample_id", "class"))
# clusters %>% filter(class1 == 2 & class2 == 1) %>% mutate(class = class3) %>% select(icgc_sample_id, class) -> a
# d_clust = rbind(d_clust, a)

for (c in classes) {
  if (cls %>% filter(class == c) %>% nrow() > 500) {
    cat("mclust on cluster: " , c, "\n")
    cls %>% filter(class == c) -> c_labels
    # data %>% filter(icgc_sample_id %in% c_labels$icgc_sample_id) -> c_data
    c_clust = subclust(c_labels$icgc_sample_id)
    d_clust = rbind(d_clust, c_clust)
  }
}

clusters %>% left_join(d_clust) %>% rename(class5 = class) -> clusters

write_tsv(clusters, "method_res/model_based/allcancer_0.001/clusters5_new.tsv")


# func --------------------------------------------------------------------

mat_name = "allcancer_0.01"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))

subclust <- function(labels) {
  data %>% filter(icgc_sample_id %in% labels) -> data1
  rownames(data1) = data1$icgc_sample_id
  mc1 = Mclust(data1[, -1])
  clust = data.frame(cbind(icgc_sample_id = data1$icgc_sample_id, class = mc1$classification), stringsAsFactors = FALSE)
  clust
}

mclust_hrch <- function(data) {
  clusters = subclust(data$icgc_sample_id)
  sample_no = length(data$icgc_sample_id)
  n = 5 # max depth
  for (d in 2:n) {
    label_d = paste0("class", d)
    cls = clusters %>% unite(class, class1:label_d, sep = ".")
    cls %>% mutate(class = sub(".NA", "", class)) -> cls
    classes = unique(cls$class)
    d_clust = data.frame(row.names = c("icgc_sample_id", "class"))
    
    for (c in classes) {
      if (cls %>% filter(class == c) %>% nrow() > 0.01 * sample_no) {
        cat("mclust on cluster: " , c, "\n")
        cls %>% filter(class == c) -> c_labels
        # data %>% filter(icgc_sample_id %in% c_labels$icgc_sample_id) -> c_data
        c_clust = subclust(c_labels$icgc_sample_id)
        d_clust = rbind(d_clust, c_clust)
      }
    }
    
    write_tsv(clusters, paste0("method_res/model_based/", mat_name, "/clusters", d, ".tsv"))
  }
}




