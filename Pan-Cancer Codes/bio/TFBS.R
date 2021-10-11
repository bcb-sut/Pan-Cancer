library(readr)
library(data.table)
library(dplyr)

tfbs_tab = fread("data/wgEncodeRegTfbsCellsV3.tab")
tfbs_bed = fread("data/wgEncodeRegTfbsClusteredV3.bed")

tfbs_bed %>% 
  select(V1:V4) %>% 
  rename(seqid = V1) %>% rename(start = V2) %>% rename(end = V3) %>% rename(TFC = V4) -> tfbs_data


files = list.files("../summarized_data/")

cons_df = data.frame(icgc_mutation_id = character(), icgc_sample_id = character(), consequence_type = character())

data %>% 
  # filter(start != end) %>% 
  mutate(pos = start) %>% select(-c(start, end)) %>% 
  inner_join(a, "seqid") %>% 
  mutate(has_tf = if_else(pos >= start & pos <= end, 1, 0)) %>%
  filter(has_tf != 0) -> tf_df


for (file in files) {
  cancer_type = sub(".tsv", "", file)
  print(cancer_type)
  data = read_tsv(paste0("../summarized_data/", file))
  
  data %>% select(icgc_mutation_id, icgc_sample_id, consequence_type) %>% distinct() -> a
  cons_df = rbind(cons_df, a)
}

write_tsv(cons_df, "data/consequence_type_df.tsv.gz")




