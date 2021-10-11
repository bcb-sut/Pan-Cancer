
library(readr)
library(plotly)
library(factoextra)
library(M3C)

mat_name = "allcancer_0.001"
data = read_tsv(file.path("matrices/sgnf_genes/", paste0(mat_name, ".tsv.gz")))
clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls = clusters %>% unite(class, class1:class2, sep = ".") %>% 
  mutate(class = ifelse(class == "2.1", paste(class, class3, sep = "."), class)) %>% select(-class3)


sgnf_0.001 = read_tsv("sgnf_genes/allcancer_0.001.tsv")
#sgnf_0.01 = read_tsv("sgnf_genes/allcancer_0.01.tsv")
cosmic_genes = read_tsv("sgnf_genes/from_dataset/census_genes.tsv")

sgnf_0.001_data <- data[ , which((names(data) %in% sgnf_0.001$gene_id) == TRUE)]


# res.pca.scaled = prcomp(data[, -1], center = T, scale. = T)

data = read_tsv("matrices/sgnf_genes/allcancer_0.001.tsv.gz")
View(data)

pca_plot <- function(data, labels, cls, plt_title) {
  # df <- data[ , which((names(data) %in% labels$gene_id) == TRUE)]
  df <- data[, -1]
  res.pca = prcomp(df)
  summary(res.pca)
  
  pca_data = as.data.frame(res.pca$x) 
  pca_data = cbind(cls, pca_data)
  
  p <- plot_ly(pca_data, x = ~PC1, y = ~PC2, color = ~factor(class), colors = "Set1") %>% 
    add_markers() %>% layout(title = plt_title)
  
  #htmlwidgets::saveWidget(p, "sgnf0.001.html")
  p
}

pca_plot(data, cosmic_genes, cls, "sgnf_0.01")

library(Rtsne)
vis <- data %>%left_join(cls, by = c("icgc_sample_id" = "icgc_sample_id")) %>% rename(label = class)

colors <- c("red", "palevioletred3","plum2", "turquoise1" , "mediumorchid1",
                                 "yellowgreen", "darkseagreen1", "darkgoldenrod1", "cyan4", 
                                 "blue4", "brown4","deeppink2",
                                 "chocolate1",  "darkslategrey","rosybrown",
                                   "springgreen4", "tan4")
vis$label = as.factor(vis$label)
numTrain <- 12270
set.seed(1)
rows <- sample(1:nrow(vis), numTrain)
train <- vis[rows,]
tsne_plot <- function(perpl=30,iterations=500,learning=200){
  set.seed(1) # for reproducibility
  tsne <- Rtsne(train[,-1], dims = 2, perplexity=perpl, verbose=TRUE, max_iter=iterations, eta=learning, check_duplicates = FALSE, pca = FALSE)
  #require(rgl)
  #plot3d(tsne$Y, col=colors[train$label])
  #legend3d("topright", legend = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"), col = colors[train$label])
  plot(tsne$Y, t='n', main = print(paste0("perplexity = ",perpl, ", max_iter = ",iterations, ", learning rate = ",learning)), xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=1, "cex.lab"=1.5, asp =1)
  text(tsne$Y, labels=".", col=colors[train$label], cex = 2.5)
}

perplexity_values <- c(40)
iteration_values <- c(1000)
learning_values <- c(1000)
theta_values <- c(0.5, 0.6)
for (p in perplexity_values){
  for( i in iteration_values){
    for( l in learning_values){
        tsne_plot(p, i, l)
    
    }
  }
}
# fviz_pca_biplot(res.pca)
# fviz_pca_contrib(res.pca.scaled)



  # layout(scene = list(xaxis = list(title = 'Weight'),
  #                     yaxis = list(title = 'Gross horsepower'),
  #                     zaxis = list(title = '1/4 mile time')))

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
# chart_link = api_create(p, filename="scatter3d-basic")
# chart_link

htmlwidgets::saveWidget(as_widget(p), "pca.html")

################################################
library(ggfortify)
library(factoextra)
feature_genes = read_tsv("sgnf_genes/gene_feature_0.001.tsv")
n_row = nrow(feature_genes)
n_col = ncol(feature_genes)
View(feature_genes[1:n_row, -1])
res.pca <- prcomp(feature_genes[1:n_row, 1:n_col], scale = TRUE)
fviz_eig(res.pca)



