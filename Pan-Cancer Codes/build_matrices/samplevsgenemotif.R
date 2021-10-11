library(readr)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SomaticSignatures)
library(GenomicRanges)

# summarized = read_tsv("../data/summarized/Bone_summarized.tsv")
dirs = list.dirs(path =
                   "../../Bayati/DeepCancer project/Data/ICGC/",
                 full.names = TRUE, recursive = FALSE)


build_genemotif_mat <- function(df) {
  if (nrow(df) == 0)
    return(data.frame())
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
    dplyr::rename(chromosome_strand = strand) %>% 
    dplyr::rename(chromosome_start = position) %>% 
    left_join(df) %>% 
    dplyr::select(icgc_sample_id, gene_affected, motif) %>% 
    drop_na() %>% 
    mutate(genemotif = paste(gene_affected, motif, sep = "_")) %>%
    group_by(icgc_sample_id, genemotif) %>% 
    summarise(count = n()) %>% 
    spread(genemotif, count, fill = 0) %>% 
    ungroup() -> sample_genemotif_mat
  # sample_motif_mat %>% 
    # mutate(total = rowSums(sample_motif_mat %>% dplyr::select(-icgc_sample_id))) -> sample_motif_mat
    
    
  return(sample_genemotif_mat)
}

sigenes = read_tsv("sig_genes/sig_genes_final.tsv")
sigenes = sigenes$gene_affected

for(dir in dirs) {
  cancer_type = sub(".*ICGC//", "", dir)
  print(cancer_type)
  data = read_tsv(paste(dir, "simple_somatic_mutation.open.tsv", sep = "/"))
  data %>% 
    filter(assembly_version == "GRCh37") %>% 
    filter(gene_affected %in% sigenes) %>% 
    dplyr::select(icgc_sample_id, chromosome, chromosome_start, chromosome_end, chromosome_strand,
                  mutation_type, mutated_from_allele, mutated_to_allele, sequencing_strategy, gene_affected) %>% 
    filter(mutation_type == "single base substitution") %>%
    dplyr::select(-mutation_type) %>% 
    # filter(chromosome_start == chromosome_end) %>%
    # filter(nchar(mutated_from_allele) == 1 & nchar(mutated_to_allele) == 1) %>% 
    # filter(mutated_from_allele != "-c" & mutated_to_allele != "-" & mutated_from_allele != "_" & mutated_to_allele != "_") %>% 
    mutate(chromosome_strand = ifelse((chromosome_strand == 1), '+', '-')) %>% 
    filter(chromosome != "MT") %>% 
    filter(mutated_from_allele != mutated_to_allele) %>%
    mutate(chromosome = paste0("chr", chromosome)) %>% 
    unique() -> df_all
  
  # df_all %>% 
  #   mutate(position = chromosome_start) %>% 
  #   mutate(strand = chromosome_strand) -> df_all
  
  df_all %>% filter(sequencing_strategy == "WGS") -> df_wgs
  df_all %>% filter(sequencing_strategy == "WXS") -> df_wxs
  
  cat("building WGS matrix for ", cancer_type, "\n")
  wgs_mat <- build_genemotif_mat(df_wgs %>% dplyr::select(-sequencing_strategy))
  write_tsv(wgs_mat,
            file.path("matrices/genemotif_sigene/", paste0(cancer_type, "_WGS_genemotif.tsv")))
  
  cat("building WXS matrix for ", cancer_type, "\n")
  wxs_mat <- build_genemotif_mat(df_wxs %>% dplyr::select(-sequencing_strategy))
  write_tsv(wxs_mat,
            file.path("matrices/genemotif_sigene/", paste0(cancer_type, "_WXS_genemotif.tsv")))
  
  # cat("merging motif matrics of ", cancer_type, "\n")
  # all_mat <- build_genemotif_mat(df_all %>% dplyr::select(-sequencing_strategy))
  # write_tsv(all_mat,
  #           file.path("matrices/genemotif_sigene/", paste0(cancer_type, "_genemotif.tsv")))
  # 
  # final_mat = plyr::rbind.fill(final_mat, sample_motif_mat)
  
}

# write_tsv(final_mat,
          # file.path("matrices/3mermotif/", "final_3mermotif.tsv"))

