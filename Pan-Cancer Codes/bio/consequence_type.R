library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

files = list.files("../summarized_data/")

cons_df = data.frame(icgc_mutation_id = character(), icgc_sample_id = character(), consequence_type = character())

for (file in files) {
  cancer_type = sub(".tsv.gz", "", file)
  print(cancer_type)
  data = read_tsv(paste0("../summarized_data/", file))
  
  data %>% select(icgc_mutation_id, icgc_sample_id, consequence_type) %>% distinct() -> a
  cons_df = rbind(cons_df, a)
}

write_tsv(cons_df, "data/consequence_type_df.tsv.gz")


# Analysis ----------------------------------------------------------------

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))
cons_df %>% 
  sample_frac(0.01) %>% 
  left_join(cls) -> cons_joined

cons_joined %>% 
  select(icgc_mutation_id, class) %>% distinct() %>% 
  group_by(class) %>% summarise(mutation_no = n()) -> mut_df # number of mutations in each class

cons_joined %>% 
  group_by(consequence_type, class) %>% 
  summarise(count = n()) %>% 
  left_join(mut_df) %>% 
  mutate(frac = count / mutation_no) -> conseq_frac

write_tsv(conseq_frac, "data/consequence_type_df.tsv")
conseq_frac = read_tsv("data/consequence_type_df.tsv")

conseq_frac$class <- factor(conseq_frac$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))

ggplot(conseq_frac %>% drop_na()) + 
  geom_bar(aes(x = consequence_type, y = frac, fill = factor(consequence_type)),
           stat = "identity", show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #theme(axis.text.x=element_blank()) + 
  # ylab("Fraction in cluster") + 
  xlab("Consequense Type") +
  facet_wrap(class ~ ., strip.position = "bottom") +
  theme(plot.tag.position = 'top',
        plot.tag = element_text()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
  labs(title = "FigureS4: Consequence type analysis",
       subtitle = "",
       caption = "",
       x = "consequence type",
       y = "fraction",
       tag = "")+
    theme(plot.title = element_text(size = 25), axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20), strip.text = element_text(size = 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
  

ggsave("bio_res/consequence_type_anaysis.pdf", width = 400, height = 350, units = "mm")

