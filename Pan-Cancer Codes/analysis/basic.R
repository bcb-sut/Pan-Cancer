library(readr)
library(tidyr)
library(mclust)

path = "matrices/sgnf_genes/"

files = list.files(path = path, full.names = FALSE, recursive = FALSE)
file = files[1]

sgnf_sample = read_tsv(file.path(path, file))
mat_name = sub(".tsv.gz", "", file)
print(mat_name)

mat = sgnf_sample
mat[is.na(mat)] <- as.numeric(0)
labels = mat$icgc_sample_id
#cancers = mat$cancer_type
# mat = mat[-c(1, 2)]
mat = mat[-1]
rownames(mat) = labels

# model based -------------------------------------------------------------

df = mat
mc = Mclust(df)

saveRDS(mc, compress = TRUE, 
        file = file.path("method_res/model_based", paste0(mat_name, ".tsv")))


# other -------------------------------------------------------------------

N = 10

fviz_nbclust(mat, kmeans, method = "wss", k.max = 15)
km.res <- kmeans(mat, N)

hc.res <- mat %>%
  scale() %>%
  eclust("hclust", k = N, graph = FALSE)

saveRDS(km.res, compress = FALSE, 
        file = file.path("method/k-means/", paste0(mat_name, ".tsv")))

saveRDS(hc.res, compress = FALSE, 
        file = file.path("method/hc/", paste0(mat_name, ".tsv")))


# dbscan ------------------------------------------------------------------

df = mc$data

dbscan::kNNdistplot(df, k = 5)

db.res <- fpc::dbscan(data = df, eps = 0.3, MinPts = 5, scale = TRUE, method = "raw")

cls = cbind(icgc_sample_id = rownames(df), class = db.res$cluster %>% as.numeric()) %>% 
  as.data.frame() %>% 
  left_join(labels_df)


hdb.res <- dbscan::hdbscan(data[, -1], minPts = 10)

hdb.res = readRDS("method_res/model_based/allcancer_0.001/hdbscan.rds")

# K-means -----------------------------------------------------------------

km.res = kmeans(mat, centers = 3)


# h -----------------------------------------------------------------------

library(apcluster)

apc.res <- apcluster(negDistMat(r = 2, data[, -1]))

