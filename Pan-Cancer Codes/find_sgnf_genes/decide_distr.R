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
gene_count = read_tsv("data/counts/gene_count.tsv.gz")
gene_count %>% dplyr::filter(cancer_type == cancer) -> df
#View(gene_count)
distr = df$count

descdist(distr, boot = 500, discrete = TRUE) 
#descdist(distr, boot = 500, discrete = FALSE)
#descdist() 
#dev.copy(descdist(distr, boot = 100, discrete = TRUE), "plots/cullen and frey.pdf")
#ggsave(file.path("plots/", paste0("cullen and frey", ".png")), width = 500, height = 200, units = "mm", limitsize = FALSE)
