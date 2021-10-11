library(readr)
library(factoextra)
library(RColorBrewer)
library(dplyr)
library(highcharter)
library(ggplot2)
library(mclust)
library(data.table)

cols <- colorRampPalette(brewer.pal(9, "Set1"))

mat_name = "allcancer_0.001"
mc = readRDS("method_res/model_based/allcancer_0.001.rds")
data = read_tsv("matrices/sgnf_genes/allcancer_0.001.tsv.gz")

fviz_mclust(mc, "BIC", palette = "jco")

fviz_mclust(mc, "classification", geom = "point", 
            pointsize = 1.5, palette = "jco")

# fviz_mclust(mc, "uncertainty", palette = "jco")

labels_df = read_tsv("data/labels.tsv")
View(labels_df)
s <-  as.data.frame( mc$classification)
setDT(s, keep.rownames = TRUE)[]
colnames(s)[1] <- "icgc_sample_id"
colnames(s)[2] <- "class"
View(s)
cls_plot(s)


# model based on each cluster ---------------------------------------------

#mc1 = subclust(1)
#saveRDS(mc1, compress = TRUE, file = file.path(paste0("method_res/model_based/", mat_name), paste0(mat_name, "_1.rds")))
#cls_plot(mc1$classification)
mc2 = subclust(2)
cls_plot(mc2$classification)


subclust <- function(clust_no) {
  mc$classification[mc$classification == clust_no] %>% as.data.frame() %>% rownames() -> c1
  data %>% filter(icgc_sample_id %in% c1) -> clust1
  rownames(clust1) = clust1$icgc_sample_id
  mc1 = Mclust(clust1[, -1])
  mc1
}

clust = data.frame(class = mc1$classification)
clust = cbind(icgc_sample_id = clust1$icgc_sample_id, clust)
# rownames(clust) = clust1$icgc_sample_id
cls_plot(clust)

cls = cbind(icgc_sample_id = clust1$icgc_sample_id, class = mc1$classification) %>% as.data.frame()
cls
cls_plot(cls)
mc2$G
summary(mc2)

# plotting ----------------------------------------------------------------

cls_plot <- function(clust) {
  # cls = data.frame(class = clust)
  # cls = data.frame(cbind(icgc_sample_id = rownames(cls), class = class))
  cls = clust %>% 
    left_join(labels_df)
  
  ngroup = length(unique(cls$cancer_type))
  plot(ggplot(cls) + geom_histogram(aes(x = factor(class), group = cancer_type, fill = cancer_type), stat = "count") + 
         scale_fill_manual(values = cols(ngroup)))
  
  cls %>% group_by(class, cancer_type) %>% summarise(count = n()) -> cls_gathered
  
  ggplot(cls_gathered, aes(y = cancer_type, x = factor(class))) +
    geom_tile(aes(fill = count)) +
    geom_text(aes(label = count)) +
         scale_fill_gradient(low = "lavender", high = "midnightblue") + theme_minimal()
  
  #View(cls_gathered)
  ggsave(file.path("bio_res", "test.pdf"), width = 250, units = "mm")
  
}


