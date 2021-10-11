library(WebGestaltR)
library(enrichR)
View(listGeneSet())
View(listIdType())
View(listReferenceSet())

resultDir = dir(path = file.path("bio_res/gene analysis/geneontology", cluster))

genes = c("") # list of interest genes
get_geneontology <- function(genes, clust) {
  dir = file.path(getwd(), paste0("bio_res/gene analysis/geneontology/", clust))
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  #enrichResult = WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
  #                         enrichDatabase = "geneontology_Biological_Process",
  #                         enrichDatabaseFile = NULL, enrichDatabaseType = NULL,
  #                         interestGene = genes, interestGeneType = "ensembl_gene_id", 
  #                         referenceSet = "genome_protein-coding",
   #                        isOutput = TRUE,
  #                         outputDirectory = dir, 
  #                         projectName = NULL)
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) {
    websiteLive <- FALSE
  }
 
  dbs <- c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
  if (websiteLive) {enriched <- enrichr(as.vector(genes), dbs)}
  #if (websiteLive) {enriched[["GO_Biological_Process_2018"]]}
  rbind(enriched[["GO_Biological_Process_2018"]], enriched[["GO_Molecular_Function_2018"]], enriched[["GO_Cellular_Component_2018"]]) %>% filter(Adjusted.P.value < 0.05) %>% select(Term, P.value, Adjusted.P.value) %>% mutate(class = clust)
}

get_pathway <- function(genes, clust) {
  dir = file.path(getwd(), paste0("bio_res/gene analysis/pathway/", clust))
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  #enrichResult = WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
   #                          enrichDatabase = "pathway_KEGG",
    #                         enrichDatabaseFile = NULL, enrichDatabaseType = NULL,
    #                         interestGene = genes, interestGeneType = "ensembl_gene_id", 
    #                        referenceSet = "genome_protein-coding",
    #                       isOutput = TRUE,
    #                         outputDirectory = dir, 
    #                         projectName = NULL)
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) {
    websiteLive <- FALSE
  }
  
  dbs <- c("KEGG_2019_Human")
  if (websiteLive) {enriched <- enrichr(as.vector(genes), dbs)}
  #if (websiteLive) {enriched[["GO_Biological_Process_2018"]]}
  View(enriched[["KEGG_2019_Human"]] )
  enriched[["KEGG_2019_Human"]] %>% filter(Adjusted.P.value < 0.05) %>% select(Term, P.value, Adjusted.P.value) %>% mutate(class = clust)
  }

get_disease <- function(genes, clust) {
  dir = file.path(getwd(), paste0("bio_res/gene analysis/disease/", clust))
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  enrichResult = WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
                             enrichDatabase = "disease_GLAD4U",
                             enrichDatabaseFile = NULL, enrichDatabaseType = NULL,
                             interestGene = genes, interestGeneType = "ensembl_gene_id", 
                             referenceSet = "genome_protein-coding",
                             isOutput = TRUE,
                             outputDirectory = dir, 
                             projectName = NULL)
  
}


# analysis ----------------------------------------------------------------

  clusts = list.dirs("bio_res/gene analysis/top genes/", full.names = FALSE)
  clusts = clusts[-1]
  
  cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
  cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
  l <- unique(cls$class)
  l[1]
  #df = read_tsv(file.path(paste0("bio_res/gene analysis/top genes"), "top100_exp_codings.tsv"))
  #View(df)
  ontology <- data.frame(term=character(), pvalue=double(),qvalue=double(), class=character())
 pathway<- data.frame(term=character(), pvalue=double(),qvalue=double(), class=character())
  
  for (i in 1:17){
    df = read_tsv(file.path(paste0("bio_res/gene analysis/top genes"), "top100_exp_codings.tsv"))
    df <- df %>% filter(clust == paste0("C", i))
    df <- df[order(df$pval), ]
    df <-head(df, 100)
    df <- df %>% left_join(hg19)
    res <-get_geneontology(df$gene_symbol , paste0("C", i))
    ontology <- rbind(res, ontology)
  
    #get_pathway(df$gene, paste0("C", i))
    #get_disease(df$gene, paste0("C", i))
  }
  View(ontology$Term %>% unique())
  write_csv(ontology, file.path("bio_res/gene analysis/geneontology", "enrichr.csv"))
  

 

  





