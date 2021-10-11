library(readr)
library(dplyr)
library(tidyr)

motifs = read_tsv("data/motifs.tsv.gz")

motifs %>% 
  group_by(icgc_sample_id, motif_3mer) %>% 
  summarise(count = n()) -> motifs_df

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls = cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp))
ngroup = unique(cls$class)
ngroup

motifs_ref <- c()
for (mutation in c('CA', 'CG', 'CT', 'TA', 'TC', 'TG')){
  for (left_base in c('A', 'C', 'G', 'T')){
    for (right_base in c('A', 'C', 'G', 'T')){
      motifs_ref <- c(motifs_ref, paste0(mutation, '-', left_base, '.', right_base))
    }
  }
}

motifs_ref_df = data.frame(motif = motifs_ref)

classes = unique(cls$class)
for (c in classes) {

  if (!dir.exists(file.path("signature", c))) {
    dir.create(file.path("signature", c))
    
    ids = cls %>% filter(class == c) %>% .$icgc_sample_id
    mat = motifs_df %>% filter(icgc_sample_id %in% ids) %>% 
      rename(motif = motif_3mer) %>%
      spread(icgc_sample_id, count, fill = 0)
    
    mat = motifs_ref_df %>% left_join(mat)
    mat[is.na(mat)] = 0
    
    write_tsv(mat[, -1], file.path("signature", c, "mat.tsv"), col_names = FALSE)
    write_tsv(data.frame(colnames(mat)[-1]), file.path("signature", c, "lables.tsv"), col_names = FALSE)
  }
}

