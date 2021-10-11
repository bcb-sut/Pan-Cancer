N_optimum <- read.table("mutational signatures analysis/output/N_optimum",  sep = " ")
clust <- c('1.1','1.2', '1.3', '1.4','1.5.1', '1.5.2','1.5.3', '1.5.4', '1.5.5', '1.5.6','1.5.7', '1.5.8', '1.6', '1.7', '1.8', '2.1', '2.2')

for (i in clust){
  df <- read_tsv(paste0("mutational signatures analysis/output/resulst_for_clust", i, "/Exposures_to_3mer_Signatures_(N=", N_optimum$V2[N_optimum$V1==i], ").tsv"))
  df <- df%>%gather(signature, exposure, 2:(N_optimum$V2[N_optimum$V1==i]+1))
  df$signature <- as.factor(df$signature)
  ggplot(df, aes(x=signature, y = exposure)) + geom_violin(trim = FALSE,  fill='#A4A4A4') +geom_boxplot(width=0.1) +  theme(axis.text.x = element_blank())+ theme_classic()
  ggsave(file.path("mutational signatures analysis/output", paste0("exposure", i, ".pdf")), width = 200, units = "mm")
}
i <- '1.1'
data <- read_tsv(paste0("mutational signatures analysis/output/resulst_for_clust", i, "/3mer_Signatures_(N=", N_optimum$V2[N_optimum$V1==i], ").tsv"))
colnames(data) <- paste0(i, colnames(data))
for (i in clust[2:17]){
  df <- read_tsv(paste0("mutational signatures analysis/output/resulst_for_clust", i, "/3mer_Signatures_(N=", N_optimum$V2[N_optimum$V1==i], ").tsv"))
  colnames(df) <- paste0(i, colnames(df))
  data <- cbind(data, df)
}
View(data)
alexandrov <- read_csv("mutational signatures analysis/output/alexandrov.csv")
alexandrov <- alexandrov[, !(colnames(alexandrov) %in% c("Type","SubType"))]
View(alexandrov)
cormat2 <-  abs(cor(data, alexandrov, method = 'pearson'))
heatmap(cormat2, scale = "row")
write.csv(data.frame(ID=row.names(cormat2),cormat2), "mutational signatures analysis/output/alex_correlation.csv", row.names=FALSE)
library(reshape2)
cormat <- abs(cor(data, method = 'pearson'))
View(colnames(cormat))
write.csv(data.frame(ID=row.names(cormat),cormat), "mutational signatures analysis/output/inter_correlation.csv", row.names=FALSE)

