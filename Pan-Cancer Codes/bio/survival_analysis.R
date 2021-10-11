{
  library(readr)
  library(survival)
  library(survminer)
  library(dplyr)
  library(ranger)
  library(ggplot2)
  library(ggfortify)
}

surv_data = read_csv("data/survival.csv")

surv_data = surv_data %>% mutate(status = ifelse(donor_vital_status == "deceased", 1, 0)) %>% 
  rename(time = donor_survival_time) %>% select(-donor_vital_status) %>% 
  arrange(class) %>% 
  # mutate(clust = paste0("C", as.numeric(as.factor(class))))
  mutate(clust = as.numeric(substring(class, 2, ))) %>% 
  mutate(time = time / 365)
surv_data$class <- factor(surv_data$class, levels = c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"))
View(surv_data$icgc_donor_id %>%unique())
ngroup = length(unique(surv_data$class))
cols <- colorRampPalette(brewer.pal(9, "Set1"))

# Kaplan Meier Analysis ---------------------------------------------------

km <- with(surv_data, Surv(time, status))
km_fit <- survfit(Surv(time, status) ~ clust, data=surv_data)
# summary(km_fit)
surv_pvalue(km_fit, combine = TRUE, get_coord = TRUE)

ggplot2::autoplot(km_fit, conf.int=F, surv.size=1, 
                  censor.size = 2 , pval = TRUE)+
scale_fill_manual(values = cols(ngroup)) + 
  labs(tag = "Figure: Survival analysis",
fill = "")+
    theme(plot.tag.position = 'top',
          plot.tag = element_text(vjust = 1, hjust = 1))
  
# facets = TRUE, ncol = 2)

ggsave(file.path("bio_res", "survival_analysis1.pdf"), width = 250, units = "mm")

ggsurv <- ggsurvplot(km_fit, data = surv_data, pval = TRUE, size= 0.5,
           legend.labs =
             c("C1", "C2", "C3", "C4","C5", "C6", "C7", "C8","C9", "C10", "C11", "C12","C13", "C14", "C15", "C16", "C17"),
           conf.int=F, surv.size=1, pval.coord = c(24, 0.8), legend.title = "Subtype", censor.size = 2.5,
            legend = "right", xlab = "Time in years")
ggsurv$plot <- ggsurv$plot + labs(
  title    = "Figure : Survival analysis",
  subtitle = "",
  tag = ""
)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
  theme(plot.tag.position = 'top',
        plot.tag = element_text(vjust = 1, hjust = 1))

ggsave(file.path("bio_res", "survival_analysis2.pdf"), width = 250, units = "mm")


