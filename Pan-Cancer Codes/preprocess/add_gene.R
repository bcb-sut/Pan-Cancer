library(readr)
library(dplyr)
library(tidyr)
library(sqldf)

files = list.files(path = "../summarized_ICGC/", 
                   full.names = FALSE, recursive = FALSE)

fantom5 = read_tsv("data/FANTOM5_gene_list.tsv")

fantom5 %>% 
  # mutate(strand = ifelse(strand == '+', 1, 0) %>% as.factor()) %>% 
  mutate(seqid = sub("chr", "", seqid)) %>% 
  rename(chromosome = seqid) -> fantom_df

for (file in files) {
  cancer_type = sub(".tsv", "", file)
  print(cancer_type)
  data = read_tsv(paste0("../summarized_ICGC/", file))
  
  dataa %>% 
    filter(mutation_type == "single base substitution") %>% 
    filter(chromosome_start == chromosome_end) %>% # omit
    mutate(position = chromosome_start) %>%
    mutate(strand = chromosome_strand) %>%  
    select(icgc_sample_id, chromosome, position, strand, mutated_from_allele, mutated_to_allele) %>% 
           # consequence_type, gene_affected, sequencing_strategy) %>% 
    # unique(chromosome = paste0("chr", chromosome)) %>% 
    mutate(strand = ifelse(strand == 1, '+', '-')) -> df
 
  # sqldf("select * from df f1 left join fantom5 f2 
  #      on (f1.chromosome == f2.seqid and f1.strand == f2.strand and f1.position >= f2.start and f1.position <= f2.end) ") -> gene_df
  # 
  gene_df_total = data.frame(row.names = "")

  # lapply(split(df, (rownames(df) %>% as.numeric()) %/% 10000), function(a) {
  #   a %>%
  #     inner_join(fantom_df) %>% # joining by chromosome and strand
  #     mutate(hasGene = ifelse(position <= end & position >= end, 1, 0)) %>%
  #     View()
  #     filter(hasGene == 1) %>%
  #     View()
  #     select(-hasGene) -> gene_df
  # })
  
  
  write_tsv(gene_df, file.path("data/gene_fantom/", paste0(cancer_type, ".tsv.gz")))
  
}

