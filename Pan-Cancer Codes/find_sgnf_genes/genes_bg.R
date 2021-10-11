{
  library(readr)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(VGAM)
  library(MASS)
  library(fitdistrplus)
  library(highcharter)
}
#using data in "data/counts/gene_count.tsv.gz" we computed the most significant genes
#for each cancer type which is obtained by using fitdistplus. the result is stored in "sgnf_genes/allcancer".

gene_count = read_tsv("data/counts/gene_count.tsv.gz")
View(gene_count)
cancers = gene_count$cancer_type %>% unique()

find_sigenes <- function(cancer) {
  print(cancer)
  gene_count %>% dplyr::filter(cancer_type == cancer) -> df
  
  # descdist(df$count, boot = 100, discrete = TRUE)
  # descdist(df$count, boot = 100, discrete = FALSE)
  
  # head and neck -> lognormal
  # esophagus -> poisson
  
  if (cancer == "Skin") {
    ## fit gamma
    fit <- fitdist(df$count, "gamma") 
    shape = coef(fit)[1]
    rate = coef(fit)[2]
    plot(ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 1000) +
           stat_function(fun = dgamma, args = list(shape = unname(shape), rate = unname(rate)), color = "blue", alpha = 0.8) +
           ggtitle(cancer))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    ggsave(file.path("plots/", paste0(cancer, ".pdf")), width = 200, units = "mm", limitsize = FALSE)
    

    
    pval = 1 - pgamma(df$count, shape = shape, rate = rate)
  }
  else {
    tb = as.data.frame(table(df$count))
    colnames(tb) = c("y", "w")
    tb$y = as.numeric(tb$y)
    tb$w = as.numeric(tb$w)
    fit <- vglm(y ~ 1, negbinomial(deviance = TRUE), data = tb,
                weights = w, crit = "coef") 
    mu = Coef(fit)[1]
    size = Coef(fit)[2]
    
    ggplot(data = df, aes(x = count)) + geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.6) + xlim(0, 1000) +
           stat_function(fun = dnbinom, args = list(size = unname(size), mu = unname(mu)), color = "blue", alpha = 0.8) +
           ggtitle(cancer)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    ggsave(file.path("plots/", paste0(cancer, ".pdf")), width = 200, units = "mm", limitsize = FALSE)
    
    
    pval = 1 - pnbinom(df$count, size = size, mu = mu)
  }
  df = cbind(df, pval)
  #View(df)
  sgnf_genes = df %>% dplyr::filter(pval < treshold) %>% dplyr::select(gene_id, cancer_type, count, pval)
  print(cat("number of sgnf in ", cancer, " : ", nrow(sgnf_genes), "of ", length(unique(df$gene_id)), "  ---------------------------------------"))
  return(sgnf_genes)
}

treshold = 0.001
sigenes_mat = data.frame(matrix(ncol = 3, nrow = 0))
colnames(sigenes_mat) <- c("gene_id", "cancer_type", "count", "pval")

for (cancer in cancers) {
  sgnf_genes = find_sigenes(cancer)
  write_tsv(sgnf_genes, file.path("plots", paste0(cancer, ".tsv")))
  sigenes_mat = rbind(sigenes_mat, sgnf_genes)
}
files = list.files(path = "data/codings/", 
                   full.names = FALSE, recursive = FALSE)

temp = data.frame(matrix(ncol = 2, nrow = 0))
colnames(temp) <- c("cancer_type", "all")

for(file in files) {
  cancer_type = sub(".tsv.gz", "", file)
  print(cancer_type)
  data = read_tsv(paste0("data/codings/", file))
  all = length(data$icgc_sample_id %>% unique())
  a <- data.frame(cancer_type, all)
  temp = rbind(temp, a)
  View(temp)
}
View(temp)
#temp <-sigenes_mat %>% group_by(cancer_type) %>%summarise(all = n())
temp2 <- sigenes_mat %>% left_join(temp, by=c("cancer_type" = "cancer_type")) %>%mutate(rate = count/all)
View(temp2)
ggplot(temp2) + geom_tile(aes(y = gene_id, x = factor(cancer_type), fill = rate)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 1)) +
  scale_fill_gradient(low = "lavender", high = "blue") + 
  theme_minimal()+
  ylab("Feature gene")+
  xlab("Cancer type") +
  labs(title = "Mutational load of feature genes in each cancer type", fill = "fraction")+
  theme(axis.text.y = element_blank(), axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14), legend.title = element_text(size = 12))
ggsave(file.path(paste0("plots/", mat_name), paste0("sgnf_genes", ".pdf")), width = 500, height = 200, units = "mm", limitsize = FALSE)


sigenes_mat %>% dplyr::select(c(gene_id, pval)) %>% group_by(gene_id) %>% summarise(pval = min(pval)) -> sgnf_genes
hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) %>% 
  dplyr::select(gene, gene_symbol) -> hg19
sgnf_genes %>% rename(gene = gene_id) %>% left_join(hg19) -> sgnf_genes
sgnf_genes %>%  dplyr::select(c(gene_symbol, pval)) -> sgnf_genes 
View(sgnf_genes)
write_tsv(sgnf_genes, file.path("sgnf_genes", "gene_feature_0.001.tsv"))