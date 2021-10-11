library(treemap)
library(data.tree)
library(readr)
library(tidyr)
library(dplyr)
library(clustree)


clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
clusters %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> clusters
df <- data.frame(clusters$class)

df <-separate(df, c("clusters.class"), into = c("c1", "c2", "c3"))
#clusters %>%  
#  group_by(class1, class2, class3, class4, class5) %>% 
#  summarise(n = n())
df$c1
df$c1 <- mapply(df$c1, FUN=as.numeric)
df$c2 <- mapply(df$c2, FUN=as.numeric)
df$c3 <- mapply(df$c3, FUN=as.numeric)


clusters$k1 = 1
clusters$k2 = df$c1
clusters$k3 = 8 * (df$c1 - 1) + df$c2
clusters$k4 = ifelse(clusters$k3 < 5, clusters$k3, ifelse(is.na(df$c3), clusters$k3 + 7, 5 + df$c3 - 1))
View(clusters)
#View(df)
#k %>% mutate(pathString = sub("/NA", "", pathString)) -> cls
#View(cls)
clustree(clusters, prefix = "k")+
  guides(edge_alpha = FALSE)+
  labs(tag = "Figure 1.b : Clustering tree",
fill = "")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 0.75, hjust = 1))
#clustree_overlay(clusters, prefix = "k", x_value = "PC1", y_value= "PC2")
ggsave(file.path("bio_res", "tree.pdf"), width = 250, height = 200 ,units = "mm")

# cls_tree <- as.Node(cls)
# 
# cls_tree$n <- function(self) sum(sapply(self$children, function(x) x$n))
# 
# print(cls_tree, "n", limit = 100)
# n
# 
# # cls_tree$Climb(position = c(1, 5, 1))$path
# 
# df = ToDataFrameTable(cls_tree, "c1", "c2", "c3", "n")
# 
# treemap(df,
#         index=c("c1", "c2", "c3"),
#         vSize="n",
#         # vColor="GNI",
#         type="value")
# 
# Count <- function(node) {
#   result <- node$n
#   if(length(result) == 0) result <- sum(sapply(node$children, Count))
#   return (result)
# }
# 
# print(cls_tree, count = Count)
# 
# SetGraphStyle(cls_tree, rankdir = "TB")
# SetEdgeStyle(cls_tree, arrowhead = "vee", color = "grey35", penwidth = 2)
# SetNodeStyle(cls_tree, style = "filled,rounded", shape = "box", fillcolor = "LightBlue", 
#              fontname = "helvetica", tooltip = GetDefaultTooltip)
# SetNodeStyle(cls_tree$n, fillcolor = "LightBlue", penwidth = "5px")
# plot(cls_tree)
# 
# cls_tree$Set(n = c(function(self) sum(sapply(self$children, 
#                                              function(child) GetAttribute(child, "n")))), 
#              filterFun = isNotLeaf)
# 
# # plot(as.dendrogram(cls_tree), center = TRUE)
# 
