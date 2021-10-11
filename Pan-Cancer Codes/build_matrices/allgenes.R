library(readr)
library(dplyr)
library(tidyr)

dirs = list.dirs(path = 
                   "../../Bayati/DeepCancer project/Data/ICGC/",
                 full.names = TRUE, recursive = FALSE)

for(dir in dirs) {
  cancer_type = sub(".*ICGC/", "", dir)
  print(cancer_type)
  data = read_tsv(paste(dir, "simple_somatic_mutation.open.tsv", sep = "/"))
  data %>% 
    # filter(consequence_type == "missense_variant") %>%
    # filter(gene_affected %in% sigenes) %>% 
    filter(mutation_type == "single base substitution") %>%
    filter(!is.na(gene_affected)) %>% 
    dplyr::select(icgc_sample_id, gene_affected, 
                  chromosome, chromosome_start, chromosome_end) %>%
    unique() %>% 
    group_by(icgc_sample_id, gene_affected) %>% 
    summarise(count = n()) %>% 
    spread(gene_affected, count, fill = 0) -> sample_gene_mat
  write_tsv(sample_gene_mat,
            file.path("matrices/allgenes_singlebase/", paste0(cancer_type, "_allgenes.tsv")))
  # final_mat = plyr::rbind.fill(final_mat, sample_gene_mat)
}

# write_tsv(final_mat,
          # file.path("res/", "final.tsv"))

