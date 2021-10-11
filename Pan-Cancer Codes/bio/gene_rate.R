library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)



######################################### pre define 2 types of data : 'data' , 'data_in_sgnf_genes'

# 'data' : a table shows number of mutation in each motif of each gene for all samples
data = read_tsv("data/sample_data.tsv")
View(data)
# feature genes : 684
sgnf_genes = read_tsv("sgnf_genes/allcancer_0.001.tsv")
View(sgnf_genes)
# a table like 'data' but only contains info about sgnf _genes in it
data_in_sgnf_gene <- data[which(data$gene_affected %in% sgnf_genes$gene_id), ]


#read clusering , the table assign each sample to a cluster
cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))
#View(cls)


hg19 = read_tsv("data/ENCODE_hg19.tsv")
hg19 %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = sub("\\.\\d+", "", gene)) -> hg19

tem <- hg19[which(hg19$Gene_class == "protein_coding"), ]
View(tem$gene_symbol %>% unique)
View(hg19)
cosmic = read_csv("data/Census_alloWed.csv")


################################################################################################

##############################fgene rate for samples have more than k mutation ########################
data_in_sgnf_gene %>% 
  select(-motif_3mer) %>% 
  left_join(cls) %>% 
  group_by(gene_affected, icgc_sample_id, class) %>%
  summarize(m_count = sum(count)) %>%
  ungroup() -> gene_joined
# k determines the threshold and gene_label can get two value : all_gen or sgnf_gene
# the function plot gene rate for all clusters
plt_gene_rate <- function(gene_label, k) {
  #if(gene_label == "All coding and lincRNA gene"){
   # gene_sample = data
  #}
  #else if(gene_label == "Significant gene"){
  #  gene_sample = data_in_sgnf_gene
  #}
  # gene_joined is a table shows number of mutations for each sample in specific genes. 
  #samples also joined with corresponding cluster here.


  #number of samples in each cluster
  gene_joined %>%
    select(-c(gene_affected, m_count)) %>% 
    group_by(class) %>%
    distinct() %>% 
    summarise(sum= n()) -> gene_sum
  
  #for each cluster and each gene shows number of samples in that cluster 
  #with at least k-mutation
  gene_joined <- gene_joined[which(gene_joined$m_count >= k ), ]
  gene_joined %>% group_by(class, gene_affected) %>% summarize(count = n()) -> gene_count
  
  #for each cluster and each gene mutation rate of gene in that cluster
  gene_count %>% left_join(gene_sum) %>% mutate(rate = count / sum) -> gene_joined
  gene_joined %>% rename(gene = gene_affected) -> gene_joined
  #protein_coding <- hg19 %>% filter(Gene_class == "lincRNA")
  sgnf_genes %>% rename(gene = gene_id) -> sgnf_genes
  merge(x=gene_joined %>% filter(class == "C1"), y = sgnf_genes, by = "gene", all.y = TRUE) -> all_gene_joined
  all_gene_joined$class = "C1"
  for (i in 2:17){
    merge(x=gene_joined %>% filter(class == paste0("C", i)), y = sgnf_genes, by = "gene", all.y = TRUE) -> temp
    temp$class = paste0("C", i)
    all_gene_joined <- rbind(all_gene_joined, temp)
  }
  all_gene_joined$rate[is.na(all_gene_joined$rate)] <- 0
  all_gene_joined$count[is.na(all_gene_joined$count)] <- 0
  all_gene_joined$sum[is.na(all_gene_joined$sum)] <- 0
  
  
  
  #View(gene_joined %>% filter(class == "C8"))
  
  ##plotting process#####
  #all_gene_joined %>% 
  #  left_join(hg19) %>%
  #  filter(Gene_class == "protein_coding") -> gene_df
    #filter(gene %in% census_genes$gene_id) 

  gene_df <- all_gene_joined

  gene_df <- gene_df %>% filter(class == "C13" | class == "C14")
  View(gene_df)
  gene_df$class <- factor(gene_df$class, levels = c("C13", "C14")) 
                                                    #"C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
  #View(gene_df$rate[gene_df $rate < 1])
  #
  #View(gene_df$rate)
  #gene_df <- gene_df[which(gene_df$class == "C14"), ]
  #write_csv(gene_df, "bio_res/gene analysis/gene rate plots/test.csv")
  #View(gene_df$rate[gene_df$rate < - 480])
  View(gene_df %>% filter(class == "C13") %>% filter(rate > 0.05))
  View(gene_df %>% filter(class == "C14") %>% filter(rate > 0.1))
  #View(gene_df %>% filter(class == "C3") %>% filter(rate > 0.78))
  #View(gene_df %>% filter(class == "C4") %>% filter(rate > 0.75))
  #View(gene_df %>% filter(class == "C6") %>% filter(rate > 0.5))
  #View(gene_df %>% filter(class == "C8") %>% filter(rate > 0.78))
  #View(gene_df %>% filter(class == "C9") %>% filter(rate > 0.78))
  #View(gene_df %>% filter(class == "C10") %>% filter(rate > 0.78))
  #View(gene_df %>% filter(class == "C11") %>% filter(rate > 0.78))
  #View(gene_df %>% filter(class == "C12") %>% filter(rate > 0.78))
  #View(gene_df %>% filter(class == "C13") %>% filter(rate > 0.95))
  #View(gene_df %>% filter(class == "C14") %>% filter(rate > 0.95))
  #View(gene_df %>% filter(class == "C15") %>% filter(rate > 0.95))
  #View(gene_df %>% filter(class == "C16") %>% filter(rate > 0.9))
  #View(gene_df %>% filter(class == "C17") %>% filter(rate > 0.95))
  
  
  
  
  p = ggplot(gene_df) +
    geom_bar(aes(x = gene, 
                  y = rate, fill = class),
             stat = "identity", show.legend = FALSE)+
    xlab("gene symbol") +
    #stat_peaks(aes(x= gene_symbol, y = rate),  color = "red", geom = "text", hjust = -0.1)+
    #stat_label_peaks(span = 41, geom = "label")+
    facet_grid(class ~ .) +
    theme(axis.text.x = element_blank()) +
    labs(tag = "", title= "Figure 2: Mutational load of feature genes")+
  ylab("fraction")+
  theme(plot.tag.position = 'top',
        plot.title = element_text(size = 30), axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30), strip.text = element_text(size = 30))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave(file.path("bio_res/gene analysis/gene rate plots", paste0(gene_label, k, ".pdf")), 
         plot = p, width = 1000, height = 600, units = "mm", dpi = 500)
  
 
}
plt_gene_rate("feature gene", 3)
#######################################################################################

############## the method computes gene rate only for samples with more than average mutation######################

plt_gene_rate_mean <- function(gene_label) {
  if(gene_label == "All coding and lincRNA gene"){
    gene_sample = data
  }
  else if(gene_label == "Significant gene"){
    gene_sample = data_in_sgnf_gene
  }
  gene_sample %>% group_by(gene_affected) %>% summarize(sum = sum(count)/12270) -> mean_mut_genes
  # gene_joined is a table shows number of mutations for each sample in specific genes. 
  #samples also joined with corresponding cluster here.
  gene_sample %>% 
    select(-motif_3mer) %>% 
    left_join(cls) %>% 
    group_by(gene_affected, icgc_sample_id, class) %>%
    summarize(m_count = sum(count)) %>%
    ungroup() -> gene_joined
  
  #number of samples in each cluster
  gene_joined %>%
    select(-c(gene_affected, m_count)) %>% 
    group_by(class) %>%
    distinct() %>% 
    summarise(sum= n()) -> gene_sum
  
  #for each cluster and each gene shows number of samples in that cluster 
  #with at least avg-mutation
  mean_mut_genes <- mean_mut_genes %>% rename(k = sum)
  gene_joined <- gene_joined %>% left_join(mean_mut_genes)
  gene_joined <- gene_joined[which(gene_joined$m_count >= gene_joined$k ), ]
  gene_joined %>% group_by(class, gene_affected) %>% summarize(count = n()) -> gene_count
  
  #for each cluster and each gene mutation rate of gene in that cluster
  gene_count %>% left_join(gene_sum) %>% mutate(rate = count / sum) -> gene_joined
  #View(gene_joined)
  
  
  ##plotting process#####
  gene_joined %>% 
    rename(gene = gene_affected) %>%
    left_join(hg19) %>% 
    filter(gene %in% census_genes$gene_id) -> gene_df
  #View(gene_df)
  
  gene_df$class <- factor(gene_df$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
  
  p = ggplot(gene_df) +
    geom_bar(aes(x = gene_symbol, 
                 y = rate, fill = class),
             stat = "identity", show.legend = FALSE)+
    xlab("gene symbol") +
    #stat_peaks(aes(x= gene_symbol, y = rate),  color = "red", geom = "text", hjust = -0.1)+
    #stat_label_peaks(span = 41, geom = "label")+
    facet_wrap(class ~ ., strip.position="bottom") +
    theme(axis.text.x = element_blank())+
    labs(tag = paste0("Figure: ", gene_label, " rate analysis with at least average", " muation"),
         ylab = "frac")+
    theme(plot.tag.position = 'top',
          plot.tag = element_text(vjust = 1, hjust = 1, size = 30), axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30), strip.text = element_text(size = 20)) 
  ggsave(file.path("bio_res/gene analysis/gene rate plots", paste0(gene_label, "average.pdf")), 
         plot = p, width = 500, height = 400, units = "mm")
}

######################################################################################################



####################### diffrential analysis plotting########################################
##########################the method plots diff analysis between each 2 subtypes ##############################
plt_dif_gene_rate <- function(gene_label, c1, c2, k) {
  if(gene_label == "All coding and lincRNA gene"){
    gene_sample = data
  }
  else if(gene_label == "Significant gene"){
    gene_sample = data_in_sgnf_gene
  }
 # View(gene_sample)
  # gene_joined is a table shows number of mutations for each sample in specific genes. 
  #samples also joined with corresponding cluster here.
  gene_sample %>% 
    select(-motif_3mer) %>% 
    left_join(cls) %>% 
    group_by(gene_affected, icgc_sample_id, class) %>%
    summarize(m_count = sum(count)) %>%
    ungroup() -> gene_joined
  
  #number of samples in each cluster
  gene_joined %>%
    select(-c(gene_affected, m_count)) %>% 
    group_by(class) %>%
    distinct() %>% 
    summarise(sum= n()) -> gene_sum
  
  #for each cluster and each gene shows number of samples in that cluster 
  #with at least k-mutation
  gene_joined <- gene_joined[which(gene_joined$m_count >= k ), ]
  gene_joined %>% group_by(class, gene_affected) %>% summarize(count = n()) -> gene_count
  
  #for each cluster and each gene mutation rate of gene in that cluster
  gene_count %>% left_join(gene_sum) %>% mutate(rate = count / sum) -> gene_joined
  View(gene_joined)
  View(sgnf_genes)
  gene_joined %>% rename(gene = gene_affected) -> gene_joined
  sgnf_genes %>% rename(gene = gene_id) -> sgnf_genes
  merge(x=gene_joined %>% filter(class == "C1"), y = sgnf_genes, by = "gene", all.y = TRUE) -> all_gene_joined
  all_gene_joined$class = "C1"
  for (i in 2:17){
    merge(x=gene_joined %>% filter(class == paste0("C", i)), y = sgnf_genes, by = "gene", all.y = TRUE) -> temp
    temp$class = paste0("C", i)
    all_gene_joined <- rbind(all_gene_joined, temp)
  }
  all_gene_joined$rate[is.na(all_gene_joined$rate)] <- 0
  all_gene_joined$count[is.na(all_gene_joined$count)] <- 0
  all_gene_joined$sum[is.na(all_gene_joined$sum)] <- 0
  
  
  
  all_gene_joined %>%
    left_join(hg19) %>%
    filter(Gene_class == "protein_coding")  -> gene_df
    #filter(gene %in% census_genes$gene_id) 
  #gene_df <- gene_joined
  #gene_df %>% rename(gene = gene_affected) -> gene_df
  gene_df %>% select(gene, class, rate) -> gene_df
  View(gene_df)

  gene_df %>% filter(class %in% c1) %>%  drop_na() %>% rename(C1 = rate) %>% select(-class) -> gene_c1
  #sgnf_genes %>% rename (gene = gene_id) %>% left_join(gene_c1, by= c("gene"= "gene")) -> gene_c1
  #gene_c1$C1[is.na(gene_c1$C1)] <- 0
  #gene_c1$class[is.na(gene_c1$class)] <- c1
  #print(length(unique(gene_c1$gene)))
  gene_df %>% filter(class %in% c2) %>% rename(C2 = rate) %>% select(-class)-> gene_c2

  #sgnf_genes  %>% rename (gene = gene_id) %>% left_join(gene_c2, by=c("gene" = "gene")) -> gene_c2
  #gene_c2$C2[is.na(gene_c2$C2)] <- 0
  #gene_c2$class[is.na(gene_c2$class)] <- c2
  #print(length(unique(gene_c2$gene)))
  #gene_c1 <- gene_c1[, !names(gene_c1) %in% "class"]
  #gene_c2 <- gene_c2[, !names(gene_c2) %in% "class"]
  
 # View(gene_c1)
#  View(gene_c2)
  
  gene_c1 %>% left_join(gene_c2) %>% mutate(dif = C1 - C2) -> diff_df
  #drops <- c("class.x","class.y")
  #diff_df <- diff_df[, !names(diff_df) %in% drops]
  diff_df  <- diff_df %>% left_join(hg19) %>% select(gene_symbol, C1, C2, dif)
  #View(diff_df)
  diff_df %>% gather("col", "val", 2:4) -> plt_df
  plt_df$col[plt_df$col == "C1"] <- c1
  plt_df$col[plt_df$col == "C2"] <- c2
  View(plt_df)
  plt_df$col = factor(plt_df$col, levels = c(c1, "dif", c2))
  p = ggplot(plt_df) +
    geom_bar(aes(x = gene_symbol, y = val), stat = "identity") +
    xlab("gene symbol") +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(axis.text.x=element_blank()) +
    facet_grid(col ~ .)+
    #stat_peaks(aes(x= gene, y = val),  color = "red", geom = "text", hjust = -0.1)+
    #stat_label_peaks(span = 41, geom = "label")+
    labs(tag = paste0("Figure : Diffrential analysis for ", c1, " & ", c2, " ", gene_label, "(k =", k , ")"),
         fill = "count",
         y = "rate")+
    ylim(-1, 1)
    theme(plot.tag.position = 'top',
          plot.tag = element_text(vjust = 1, hjust = 1))
  
  
  ggsave(file.path("bio_res/gene analysis/differential analysis", paste0(gene_label, k,  " ", c1, " & ", c2, ".pdf")), 
         plot = p)
          #width = 200, height = 500, units = "mm")
}
#############################################################################################################

plt_dif_gene_rate("Significant gene", "C16", "C17", 1)

######################## compute dif value between two clusters #############################3
dif_vs_compute <-function(gene_label, c1, c2, k) {
  if(gene_label == "All coding and lincRNA gene"){
    gene_sample = data
  }
  else if(gene_label == "Significant gene"){
    gene_sample = data_in_sgnf_gene
  }
  # View(gene_sample)
  # gene_joined is a table shows number of mutations for each sample in specific genes. 
  #samples also joined with corresponding cluster here.
  gene_sample %>% 
    select(-motif_3mer) %>% 
    left_join(cls) %>% 
    group_by(gene_affected, icgc_sample_id, class) %>%
    summarize(m_count = sum(count)) %>%
    ungroup() -> gene_joined
  
  #number of samples in each cluster
  gene_joined %>%
    select(-c(gene_affected, m_count)) %>% 
    group_by(class) %>%
    distinct() %>% 
    summarise(sum= n()) -> gene_sum
  
  #for each cluster and each gene shows number of samples in that cluster 
  #with at least k-mutation
  gene_joined <- gene_joined[which(gene_joined$m_count >= k ), ]
  gene_joined %>% group_by(class, gene_affected) %>% summarize(count = n()) -> gene_count
  
  #for each cluster and each gene mutation rate of gene in that cluster
  gene_count %>% left_join(gene_sum) %>% mutate(rate = count / sum) -> gene_joined
  
  
  gene_df <- gene_joined
  gene_df %>% rename(gene = gene_affected) -> gene_df
  gene_df %>% select(gene, class, rate) -> gene_df
  # View(gene_df)
  
  gene_df %>% filter(class %in% c1) %>%  drop_na() %>% rename(C1 = rate) %>% select(-class) -> gene_c1
  sgnf_genes %>% rename (gene = gene_id) %>% left_join(gene_c1, by= c("gene"= "gene")) -> gene_c1
  gene_c1$C1[is.na(gene_c1$C1)] <- 0
  gene_c1$class[is.na(gene_c1$class)] <- c1
  #print(length(unique(gene_c1$gene)))
  gene_df %>% filter(class %in% c2) %>% rename(C2 = rate) %>% select(-class)-> gene_c2 
  sgnf_genes  %>% rename (gene = gene_id) %>% left_join(gene_c2, by=c("gene" = "gene")) -> gene_c2
  gene_c2$C2[is.na(gene_c2$C2)] <- 0
  gene_c2$class[is.na(gene_c2$class)] <- c2
  #print(length(unique(gene_c2$gene)))
  gene_c1 <- gene_c1[, !names(gene_c1) %in% "class"]
  gene_c2 <- gene_c2[, !names(gene_c2) %in% "class"]
  
  # View(gene_c1)
  #  View(gene_c2)
  
  gene_c1 %>% left_join(gene_c2) %>% mutate(dif = C1 - C2) -> diff_df
  
  #drops <- c("class.x","class.y")
  #diff_df <- diff_df[, !names(diff_df) %in% drops]
  #View(diff_df)"
  diff_df <- diff_df[, !names(diff_df) %in% c("C1", "C2")]
  diff_df
}
  #################################################################################

######################### table for dif value #####################################
dif_table_func <- function(gene_label, k, c){
  table_c  <- dif_vs_compute(gene_label = gene_label, c1= c, c2 = "C1", k = k)
  names(table_c)[names(table_c) == "dif"] <- "C1"
  for(i in 2:17){
    temp <- dif_vs_compute(gene_label = gene_label, c1= c, c2 = paste0("C", i) , k = k)
    names(temp)[names(temp) == "dif"] <- paste0("C", i)
    table_c <- table_c %>% left_join(temp, by = c("gene" = "gene"))
    View(table_c)
  }
  write_csv2(table_c, file.path("bio_res/gene analysis/differential analysis/", paste0(c, ".csv")))
}
#########################################################################################


cosmic %>% 
  select(Synonyms) %>% 
  rename(gene_id = Synonyms) %>% 
  separate_rows(gene_id, sep = ",") %>% 
  filter(gene_id %in% hg19$gene) -> census_genes
View(cosmic)

#plt_gene_rate("Significant gene", 1)
#plt_gene_rate("Significant gene", 2)
#plt_gene_rate("Significant gene", 5)
#plt_gene_rate("All coding and lincRNA gene", 1)
#plt_gene_rate("All coding and lincRNA gene", 2)
#plt_gene_rate("All coding and lincRNA gene", 5)
plt_dif_gene_rate("Significant gene", "C1", "C2", 5)
for(i in  1:17){
  for (j in (i+1):17){
    #dif_func("All coding and lincRNA gene", paste0("C", i), paste0("C", j), 1)
    #dif_func("All coding and lincRNA gene", paste0("C", i), paste0("C", j), 2)
    #dif_func("All coding and lincRNA gene", paste0("C", i), paste0("C", j), 5)
    plt_dif_gene_rate("Significant gene", paste0("C", i), paste0("C", j), 1)
    plt_dif_gene_rate("Significant gene", paste0("C", i), paste0("C", j), 2)
    plt_dif_gene_rate("Significant gene", paste0("C", i), paste0("C", j), 5)
  }
}

for (i in 1:17){
  dif_table_func("Significant gene", 1, paste0("C", i))
}

#plt_gene_rate_mean("Significant gene")
#plt_gene_rate_mean("All coding and lincRNA gene")

