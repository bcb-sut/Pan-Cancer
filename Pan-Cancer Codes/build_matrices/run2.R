library(readr)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SomaticSignatures)
library(GenomicRanges)

dirs = list.dirs(path = 
                   "../../Bayati/DeepCancer project/Data/ICGC/",
                 full.names = TRUE, recursive = FALSE)

final_mat = data.frame()

for(dir in dirs) {
  data = read_tsv(paste(dir, "simple_somatic_mutation.open.tsv", sep = "/"))
  cancer_type = sub(".*ICGC//", "", dir)
  print(cancer_type)
  data %>% 
    filter(consequence_type == "missense_variant") %>% 
    dplyr::select(icgc_sample_id, chromosome, chromosome_start, chromosome_end, chromosome_strand,
                  mutated_from_allele, mutated_to_allele) %>% 
    filter(nchar(mutated_from_allele) == 1 & nchar(mutated_to_allele) == 1) %>% 
    filter(mutated_from_allele != "-" & mutated_to_allele != "-" & mutated_from_allele != "_" & mutated_to_allele != "_") %>% 
    filter(chromosome != "X") %>%
    mutate(chromosome_strand = ifelse((chromosome_strand == 1), '+', '-')) %>% 
    mutate(chromosome = paste0("chr", chromosome)) -> df
  
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns = T,
                                 seqnames.field = 'chromosome',
                                 start.field ='chromosome_start',
                                 end.field = 'chromosome_end',
                                 strand.field = 'chromosome_strand')
  
  idx <- match(c('mutated_from_allele', 'mutated_to_allele', 'icgc_sample_id'), names(mcols(gr)))
  
  mcols(gr) <- cbind(mcols(gr)[idx], mcols(gr)[-idx])
  
  vr <- makeVRangesFromGRanges(gr, ref.field='mutated_from_allele',
                               alt.field='mutated_to_allele',
                               sampleNames.field = 'icgc_sample_id',
                               keep.extra.columns = T)
  
  vr <- mutationContext(vr, Hsapiens)
  
  
  variations <- data.frame(icgc_sample_id = mcols(gr)$icgc_sample_id,
                           chromosome = as.character(seqnames(vr)),
                           position = start(vr),
                           strand = as.character(strand(gr)),
                           motif = paste0(as.character(mcols(vr)$alteration), '-', as.character(mcols(vr)$context)))
  variations %>% 
    unique() %>% 
    group_by(icgc_sample_id, motif) %>% 
    summarise(count = n()) %>% 
    spread(motif, count, fill = 0) -> sample_motif_mat
  write_tsv(sample_motif_mat,
            file.path("res/", paste0(cancer_type, "_motif.tsv")))
  
  final_mat = plyr::rbind.fill(final_mat, sample_motif_mat)
  
}

write_tsv(final_mat,
          file.path("res/", "final_motif.tsv"))

