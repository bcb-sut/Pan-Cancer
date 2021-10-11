library(readr)
library(dplyr)

specimen = read_tsv("../specimen.tsv.gz")
somatic = read_tsv("../simple_somatic_mutation.open.tsv.gz")

normal_labels = specimen %>% filter(startsWith(specimen_type, "Normal")) %>% .$icgc_specimen_id
tumor_labels = specimen %>% filter(startsWith(specimen_type, "Primary tumour")) %>% .$icgc_specimen_id

files = list.files(path = "../summarized_ICGC/", 
                   full.names = FALSE, recursive = FALSE)

for (file in files) {
  cancer_type = sub(".tsv", "", file)
  print(cancer_type)
  data = read_tsv(paste0("../summarized_ICGC/", file))
  data %>% 
    filter(mutation_type == "single base substitution") %>% 
    filter(chromosome_start == chromosome_end) %>% # omit
    mutate(position = chromosome_start) %>%
    mutate(strand = chromosome_strand) %>%  
    select(icgc_sample_id, icgc_specimen_id, chromosome, position, strand, mutated_from_allele, mutated_to_allele) %>% 
    # consequence_type, gene_affected, sequencing_strategy) %>% 
    # unique(chromosome = paste0("chr", chromosome)) %>% 
    mutate(strand = ifelse(strand == 1, '+', '-')) -> df
  
  df %>% filter(icgc_specimen_id %in% tumor_labels) -> tumors_df
  df %>% filter(icgc_specimen_id %in% normal_labels) -> normal_df
  
  gene_df_total = data.frame(row.names = "")
  
  write_tsv(gene_df, file.path("data/gene_fantom/", paste0(cancer_type, ".tsv.gz")))
  
}

