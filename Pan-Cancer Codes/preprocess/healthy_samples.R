library(readr)
library(vcfR)
library(sqldf)
library(dplyr)
library(tidyr)

# ch22 = read.vcfR(limit = 1e10, "../../Rezaei/input/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
        
file = "../../Rezaei/input/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

tmp.vcf = readLines(file)
data = read.table(gzfile(file))

tmp.vcf = tmp.vcf[-(grep("#CHROM",tmp.vcf) + 1) : -(length(tmp.vcf))]
vcf.names = unlist(strsplit(tmp.vcf[length(tmp.vcf)], "\t"))
names(data) = vcf.names

# write_tsv(tmp.vcf.data, "data/healthy_samples/chr22.tsv.gz")

data = read_tsv("data/healthy_samples/chr22.tsv.gz")
data = tmp.vcf.data
colnames(data)[1] = "CHROM"
chrom_no = 22

fantom5 = read_tsv("data/FANTOM5_gene_list.tsv")

fantom5 %>%
  filter(gene_biotype == "protein_coding") %>% 
  mutate(seqid = sub("chr", "", seqid)) %>% 
  filter(seqid == chrom_no) %>% 
  mutate(seqid = as.numeric(seqid)) %>% 
  select(-gene_biotype) %>% 
  rename(CHROM = seqid)-> fantom_df

data %>%
  dplyr::select(CHROM, POS) -> df

# sqldf("select * from df f1 left join fantom_df f2 
#        on (f1.CHROM == f2.seqid and f1.POS >= f2.start and f1.POS <= f2.end)") %>% 
#   # select(CHROM, POS, gene_id) %>% 
#   drop_na() -> gene_df

gene_df_total = data.frame(row.names = "POS", "seqid", "strand", "gene_id", "REF", "ALT", "sample_id", "mutation")

lapply(split(df, (rownames(df) %>% as.numeric()) %/% 100), function(a) {
  a %>%
    inner_join(fantom_df) %>% # joining by chromosome
    mutate(hasGene = ifelse(POS <= end & POS >= end, 1, 0)) %>%
    filter(hasGene == 1) %>%
    select(-hasGene) -> gene_df
  gene_df %>%
    left_join(data) %>% 
    select(-c(CHROM, start, end)) %>% 
    select(-c(ID, QUAL, FILTER, INFO, FORMAT)) %>% 
    gather(HG00096:NA21144, key = "sample_id", value = "mutation") %>% 
    filter(mutation != "0|0") -> gene_sample
  gene_df_total = rbind(gene_df_total, gene_sample)
 
})

  
write_tsv(gene_df_total, file.path("../healthy_samples/", paste0("chr", chrom_no, ".tsv")))

