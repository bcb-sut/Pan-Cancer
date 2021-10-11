library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters3.tsv")
cls = clusters %>% unite(class, class1:class2, sep = ".") %>% 
  mutate(class_sys = ifelse(class == "2.1", paste(class, class3, sep = "."), class)) %>% 
  mutate(class_exp = ifelse(class == "1.5", paste(class, class3, sep = "."), class)) %>% 
  select(-c(class, class3))

# write_tsv(cls, "method_res/model_based/allcancer_0.001/clusters.tsv")

samplevsmotif = read_csv("data/3mermotif.csv")

sig_data = read_csv("data/signatures_probabilities.csv") %>% select(-c(X34:X40))
sig_data %>% mutate(motif = gsub('[[:punct:] ]+' , '', `Somatic Mutation Type`)) %>% 
  mutate(from_to = substr(motif, 2, 3)) %>%
  mutate(context = paste0(substr(motif, 1, 1), substr(motif, 4, 4))) -> sig

samplevsmotif %>% left_join(cls) %>% 
  gather(`CA-A.A`:`CA-A.N`, key = "motif", value = "rate") %>% 
  drop_na() %>% 
  group_by(class_sys, motif) %>% 
  dplyr::mutate(rate = as.numeric(rate)) %>% 
  dplyr::summarise(rate = sum(rate)) %>% 
  mutate(motif = gsub('[[:punct:] ]+' , '', motif)) -> df_samplemotif

df_samplemotif %>% 
  group_by(class_sys) %>% 
  mutate(rate = (rate / sum(rate))) %>%
  mutate(from_to = substr(motif, 1, 2)) %>%
  mutate(context = substr(motif, 3, 4)) %>% 
  select(-motif) %>%
  spread(class_sys, rate, fill = 0) -> sig2

sig %>% left_join(sig2) %>% 
  select(-c(motif, `Substitution Type`, Trinucleotide, `Somatic Mutation Type`, from_to, context)) -> a



a = as.matrix(a)
t<-matrix(0,nrow = 16 ,ncol = 30)
a1<-a[,1:30]
a2<-a[,31:46]
for (i in 1:16) {
  for (j in 1:30) {
      t[i,j] <- lsa::cosine((a2[ ,i]), (a1[ ,j]))
  }
}

# vec1 = a1[,1]
# vec2 = a2[,1]
# lsa::cosine(vec1,vec2)
library(RColorBrewer)
# coul = colorRampPalette(brewer.pal(8, "skyblue"))
rownames(t) = colnames(a)[31:46]
heatmap(t, Colv = NA, Rowv = NA, xlab = "signatures", ylab = "Clusters", 
        col= colorRampPalette(brewer.pal(8, "Blues"))(25))


