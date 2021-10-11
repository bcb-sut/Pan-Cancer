library(readr)
library(dplyr)
library(tidyr)

labels = read_tsv("res/sample_cancer.tsv")
mat = read_tsv("matrices/allgenes.tsv")

labels %>% 
  left_join(mat) %>% 
  drop_na() -> mat_label

write_tsv(mat_label, "matrices/allgenes.tsv.gz")

