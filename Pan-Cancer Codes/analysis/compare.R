library(fpc)

distance = read_tsv("matrices/sgnf_genes/allcancer_0.001_dist.tsv.gz")
clusters = read_tsv("method_res/model_based/allcancer_0.001/clusters3.tsv")

sys.res = clusters %>% unite(class, class1:class2, sep = ".") %>% 
  mutate(class = ifelse(class == "2.1", paste(class, class3, sep = "."), class)) %>% select(-class3) %>% 
  mutate(sys_cluster = as.numeric(as.factor(class)))

emp.res = clusters %>% unite(class, class1:class2, sep = ".") %>% 
  mutate(class = ifelse(class == "1.5", paste(class, class3, sep = "."), class)) %>% select(-class3) %>% 
  mutate(emp_cluster = as.numeric(as.factor(class)))

# comparing 2 cluster solutions
stt = cluster.stats(distance, sys.res$cluster, sys.res$cluster)
saveRDS(stt, "matrices/sgnf_genes/stt_0.001.rds")

df = left_join(sys.res %>% select(icgc_sample_id, sys_cluster), 
               emp.res %>% select(icgc_sample_id, emp_cluster), by = "icgc_sample_id") %>% 
  group_by(sys_cluster, emp_cluster) %>% 
  summarise(count = n())

ggplot(df) + geom_tile(aes(y = factor(sys_cluster), x = factor(emp_cluster), fill = count)) +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  # theme_minimal() +
  xlab("Systematic Cluster No. ") + ylab("Emperical Cluster No. ")
