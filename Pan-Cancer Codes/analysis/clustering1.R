library(readr)
library(factoextra)

allgenes_singlebase = read_tsv("matrices/allgenes.tsv.gz")

mat = allgenes_singlebase

mat[is.na(mat)] <- 0
labels = mat$icgc_sample_id
mat = mat[-1]
rownames(mat) = labels

dist = dist(mat, method = "euclidean")

clus_centroid = hclust(dist, method = "centroid")
clus = clus_centroid
plot(clus)
rect.hclust(clus, 9)

kcl = kmeans(mat, centers = 9)
fviz_cluster(kcl, mat)

plot(mat, col = kcl$cluster, pch = 20, cex = 2)


# ----

pca = prcomp(mat)
fviz_pca_biplot(pca, habillage = as.factor(kcl$cluster))




