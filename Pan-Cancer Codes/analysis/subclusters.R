{
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
}

labels_df = read_tsv("data/labels.tsv")

# number of patients in each cancer type

labels_df %>% 
  group_by(cancer_type) %>% 
  summarise(n = n()) -> cancer_count

cls_plot <- function(clust, plot_name) {
  cls = clust %>% 
    left_join(labels_df) %>% as.data.frame()
  
  ngroup = length(unique(cls$cancer_type))
  p1 <- ggplot(cls) + geom_histogram(aes(x = factor(class), fill = cancer_type), stat = "count") + 
    scale_fill_manual(values = cols(ngroup)) + 
    xlab("class") + ggtitle(plot_name)
  plot(p1)
  ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_bar.png")))
  cls %>% group_by(class, cancer_type) %>% summarise(count = n()) %>% left_join(cancer_count) -> cls_gathered
  
  p2 <- ggplot(cls_gathered) + geom_tile(aes(y = cancer_type, x = factor(class), fill = (count / n))) +
    scale_fill_gradient(low = "white", high = "steelblue") + ggtitle(plot_name) + xlab("class")
  plot(p2)
  ggsave(file.path(paste0("../plots/", mat_name), paste0(plot_name, "_heatmap.png")))
}
plot_subclusts <- function(data, mc, mcarray, plot_name, clust_nos) {
  clust = data.frame(cbind(icgc_sample_id = data$icgc_sample_id, class = mc$classification), stringsAsFactors = FALSE)
  n = length(unique(clust$class))
  for (i in clust_nos) {
    clust %>% filter(class == i) %>%  .$icgc_sample_id -> c1
    clust1 = cbind(icgc_sample_id = c1, class = mcarray[[i]]$classification) %>% as.data.frame()
    cls_plot(clust1, paste(plot_name, i, sep = "."))
    cat("cluster no. ", i, "\t number of subclusters: ", mcarray[[i]]$G, "\n")
  }
}

mat_name = "allcancer_0.001"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))

files = list.files(file.path("method_res/model_based", mat_name), full.names = TRUE)
mc = readRDS("method_res/model_based/allcancer_0.001.rds")
mc1 = readRDS(files[1])
mc2 = readRDS(files[2])
mcarray <- vector("list", 2)
mcarray[[1]] = readRDS(files[1])
mcarray[[2]] = readRDS(files[2])



clust = data.frame(cbind(icgc_sample_id = data$icgc_sample_id, class = mc$classification), stringsAsFactors = FALSE)
plot_name =  paste(mat_name, " ")
cls_plot(clust, plot_name)
plot_subclusts(data, mc, mcarray, paste0(plot_name, "c"))

subclust <- function(mc, data, clust_nos) {
  mcarray = vector("list", mc$G)
  clust = data.frame(cbind(icgc_sample_id = data$icgc_sample_id, class = mc$classification), stringsAsFactors = FALSE)
  for (cluster_no in clust_nos) {
    clust %>% filter(class == cluster_no) %>%  .$icgc_sample_id -> c1
    data %>% filter(icgc_sample_id %in% c1) -> clust1
    rownames(clust1) = clust1$icgc_sample_id
    mc1 = Mclust(clust1[, -1])
    mcarray[[cluster_no]] = mc1
  }
  mcarray
}

clust = data.frame(cbind(icgc_sample_id = data$icgc_sample_id, class = mc$classification), stringsAsFactors = FALSE)
clust %>% filter(class == 1) %>% .$icgc_sample_id -> c1
data %>% filter(icgc_sample_id %in% c1) -> clust1
clust %>% filter(class == 2) %>% .$icgc_sample_id -> c2
data %>% filter(icgc_sample_id %in% c2) -> clust2

mcarr1 = subclust(mc1, clust1, c(5))
plot_subclusts(clust1, mc1, mcarr1, paste0(plot_name, "c 1"), c(5))
mcarr2 = subclust(mc2, clust2,  c(1, 2))
