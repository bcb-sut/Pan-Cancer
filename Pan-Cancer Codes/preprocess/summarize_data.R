library(readxl)
library(data.table)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(readr)

dirs = list.dirs(path = 
                   "../../ICGC_data_2017/",
                 # "../../../data/dashti/",
                 full.names = TRUE, recursive = FALSE)

for (dir in dirs) {
  cancer_type = sub(".*2017//", "", dir)
  data = read_tsv(paste(dir, "simple_somatic_mutation.open.tsv.gz", sep = "/"))
  
  data %>% 
    filter(assembly_version == "GRCh37") %>% 
    filter(mutation_type == "single base substitution") %>%
    filter(chromosome != "MT") %>% 
    mutate(seqid = paste0("chr", chromosome)) %>% 
    mutate(start = chromosome_start) %>% 
    mutate(end = chromosome_end) %>% 
    # filter(start == end) %>% 
    mutate(strand = chromosome_strand) %>% 
    select(icgc_mutation_id, icgc_sample_id, seqid, start, end, strand, gene_affected,
           mutated_from_allele, mutated_to_allele, sequencing_strategy, consequence_type) %>% 
    distinct() -> mutations_loci_samples
  
  write_tsv(mutations_loci_samples, 
            file.path("../summarized_data/", paste0(cancer_type, ".tsv.gz")))
}

