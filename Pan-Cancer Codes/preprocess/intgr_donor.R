library(readr)
library(dplyr)

dirs = list.dirs(path = 
                   "../../ICGC_data_2017/",
                 full.names = TRUE, recursive = FALSE)

# sample_df = data.frame(icgc_sample_id = character(), icgc_donor_id = character())

donor_df = data.frame(icgc_sample_id = character(), icgc_donor_id = character(), project_code = character(),
                      donor_sex = character(), donor_vital_status = character(), donor_survival_time = integer(),
                      donor_age_at_diagnosis = integer(), donor_age_at_enrollment = numeric())

for (dir in dirs) {
  cancer_type = sub(".*2017//", "", dir)
  sample = read_tsv(paste(dir, "sample.tsv.gz", sep = "/"))
  donor = read_tsv(paste(dir, "donor.tsv.gz", sep = "/"))
  
  d = donor %>% select(icgc_donor_id, project_code, donor_sex, donor_vital_status, donor_survival_time,
                   donor_age_at_diagnosis, donor_age_at_enrollment)
  s = sample %>% select(icgc_sample_id, icgc_donor_id)
  
  donor_df = rbind(donor_df, s %>% left_join(d))
  
}

write_tsv(donor_df, "data/donor_meta.tsv.gz")

