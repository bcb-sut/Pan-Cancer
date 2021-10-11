library(readr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(tidyverse)

# cols <- colorRampPalette(brewer.pal(9, "Set1"))

meta_data = read_tsv("data/donor_meta.tsv")

cls = read_tsv("method_res/model_based/allcancer_0.001/clusters.tsv")
cls %>% mutate(class = class_exp) %>% select(-c(class_sys, class_exp)) -> cls
cls$class = paste0("C", as.character(as.numeric(as.factor(cls$class))))

donor_data = left_join(cls, 
                       select(meta_data, 
                              icgc_sample_id, icgc_donor_id, project_code, donor_sex, donor_vital_status, donor_survival_time) %>% 
                         distinct())

cls %>% group_by(class) %>% summarize(inclass = n()) -> n_class

donor_data %>% 
  select(icgc_sample_id, class, project_code) %>% 
  distinct() %>% 
  mutate(country = sub(".*-", "", project_code)) %>% 
  group_by(class, country) %>% 
  summarize(count = n()) %>% 
  left_join(n_class, "class") %>% 
  mutate(freq = count / inclass) -> project_df_country



project_df_country$class <- factor(project_df_country$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
# ggplot(project_df, aes(x = country, y = class, fill = freq)) + 
#   geom_tile(color = "grey") + 
#   # scale_fill_gradient2(low="red", high="blue") +
#   # theme_dark() +
#   # xlab('Project Code') +
#   xlab('Region') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + 
#   coord_fixed()
perm <- expand.grid(unique(project_df_country$class), unique(project_df_country$country)) %>%rename (class = Var1) %>%rename(country = Var2)
View(perm)
perm %>% left_join(project_df_country)-> project_df_country
project_df_country$freq[is.na(project_df_country$freq)] <- 0 
project_df_country$freq[project_df_country$freq < 0.01] <- 0
countrymat <- as.matrix(project_df_country)
countrymat <- sweep(countrymat, 1, rowSums(countrymat))
project_df_country <- as.data.frame(countrymat)


ggplot(project_df_country, aes(x = country, y = class)) + 
  geom_tile(aes(fill = freq), size = 8) +
  geom_text(aes(label = round(freq, 2)), size = 3, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) +
  xlab('Region')+
  labs(title = "Figure 4.b : Region distribution",
       subtitle = "",
       caption = "",
       x = "Region", y = "Subtype",
       tag = "",
       fill = "fraction")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))

ggsave("bio_res/region_analysis.pdf", width = 250, units = "mm")

donor_data %>% 
  select(icgc_sample_id, class, project_code) %>% 
  distinct() %>% 
  group_by(class, project_code) %>% 
  summarize(count = n()) %>% 
  left_join(n_class, "class") %>% 
  mutate(freq = count / inclass) -> project_df
project_df$class <- factor(project_df$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(project_df, aes(x = project_code, y = class)) + 
  geom_tile(aes(fill = freq), color = "grey", size = 4) +
  geom_text(aes(label = round(freq, 2)), size = 3) +
  scale_fill_gradient(low = "lavender", high = "midnightblue") +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) +
  xlab('Project code')+
  labs(title = "",
       subtitle = "",
       caption = "",
       x = "Project code", y = "Subtype",
       tag = "Supplementary Figure : Project code distribution of samples in each subtypes",
       fill = "fraction")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggsave("bio_res/projectcode_analysis.pdf", width = 1000, units = "mm")



# gender ------------------------------------------------------------------

donor_data %>% 
  select(icgc_sample_id, class, donor_sex) %>% 
  distinct() %>% 
  group_by(class, donor_sex) %>% 
  summarize(count = n()) %>% 
  left_join(n_class, "class") %>% 
  mutate(freq = count / inclass) -> gender_df
gender_df <- gender_df %>% drop_na()

gender_sum_df <- gender_df %>% group_by(class) %>% summarize(sum = sum(freq))
gender_df <- gender_sum_df %>% left_join(gender_df, "class") %>%mutate(freq = freq/sum)

View(gender_df)
gender_df$class <- factor(gender_df$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
ggplot(gender_df %>% drop_na(), aes(x = donor_sex, y = class)) + 
  geom_tile(aes(fill = freq), color = "grey") +
  geom_text(aes(label = round(freq, 2)), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) +
  xlab('Gender')+
  labs(title = "Figure 4.a : Gender distribution",
       subtitle = "",
       caption = "",
       x = "Gender", y = "Subtype",
       tag = "",
       fill = "fraction")+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))
  

ggsave("bio_res/gender_analysis.pdf", width = 180, units = "mm")

# survival  ---------------------------------------------------------------

write_csv(donor_data %>% 
            select(icgc_donor_id, class, donor_vital_status, donor_survival_time) %>% drop_na(), "data/survival.csv")



