#!/usr/bin/env Rscript

require(bigrquery)
require(tidyverse)
require(dplyr)
require(ggplot2)
require(glue)

source("/mnt/arcus/lab/users/kafadare/data_functions.R")

#ICD to phecode conversion tables downloaded from https://phewascatalog.org/ on Sep 6 2024
icd10tophecode <- read.csv('/mnt/arcus/lab/users/kafadare/phecode_icd10.csv')
icd9tophecode <- read.csv('/mnt/arcus/lab/users/kafadare/phecode_icd9_rolled.csv') 

phecode_defs <- read.csv('/mnt/arcus/lab/users/kafadare/filter-by-dx/dx-filters/excuded_phecodes2.csv')
phecodes_exclusion <- read.csv('/mnt/arcus/lab/users/kafadare/filter-by-dx/dx-filters/phecodes_with_exclusion_TS_and_AAB_19April2024.csv')

excluded_phecodes <- phecodes_exclusion[which(phecodes_exclusion$exclude_or_include_AAB_TS == 'exclude'),]
icd10Dx_toexclude <- icd10tophecode[which(icd10tophecode$PheCode %in% excluded_phecodes$phecode),]
icd9Dx_toexclude <- icd9tophecode[which(icd9tophecode$PheCode %in% excluded_phecodes$phecode),]

n <- as.numeric(Sys.getenv (x = "SPLIT_ID"))
print(paste0("SPLIT_ID is: ", n))

#read in the arguments passed on from the shell script
args <- commandArgs(trailingOnly = TRUE)

# Is argument length > 0?
if (length(args) > 0) {
  chunks <- as.numeric(args[1])
}else{
  n = 1
  chunks = 1
  print(paste0("SPLIT_ID is: ", n))
}


print(paste0("NUMBER OF CHUNKS is: ", chunks))

bq_auth()
proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]

#Look at neuropsych diagnoses in the included COHORT
phecodes_of_interest <- c(315.1, 315.2, 313.1, 313.3, 300.1, 300.3, 296, 296.1, 296.22, 295, 295.1, 295.2, 295.3) 

#load all the saved include-dx dataframes
path <- "/mnt/arcus/lab/users/kafadare/GAB_cohort/included/"
included_df_files <- list.files(path = path, pattern = ".*\\.rds$", full.names = TRUE)
print(paste0("Reading in files from: ", path))
included_df_list <- lapply(list.files(path = path, pattern = ".*\\.rds$", full.names = TRUE), readRDS)
print(paste0("Files read, dimensions of the list: ", length(included_df_list)))

columns_of_interest <- c("dx_id", "epic_dx_id", "age_in_days", "proc_ord_year", 
                      "proc_name", "proc_ord_id", "proc_ord_datetime","proc_id",
                      "pat_id", "sex", "pat_status_code", "pat_status_name",
                      "dob", "dod", "race", "ethnicity", "religion", "marital_status",
                      "gender_identity_code", "gender_identity_name", "language",
                      "birth_length_in", "birth_length_cm","birth_weight_oz", "birth_weight_kg",
                      "gestational_age", "gestational_age_num", "country", "state", "state_abbr",
                      "partial_zip", "opted_out_ind", "seen_last_two_yrs_ind", "deceased_ind",
                      "dob_year", "dod_age", "dod_year", "protected_ind", "extracted_date",
                      "procedure_datetime", "dx_mri_deltaT", "phecode", "description", "group",
                      "exclude_flag_AABTS", "age_at_scan_months", "icd10_list", "icd9_list")

columns_to_merge <- c("proc_ord_id",
                         "pat_id", "sex",
                         "dob", "dod", "race", "ethnicity", "religion",
                         "gender_identity_code", "gender_identity_name", "language",
                         "birth_length_in", "birth_length_cm","birth_weight_oz", "birth_weight_kg",
                         "gestational_age", "gestational_age_num", "country", "state", "state_abbr",
                         "partial_zip", "opted_out_ind", "seen_last_two_yrs_ind", "deceased_ind",
                         "dob_year", "dod_age", "dod_year", "protected_ind",
                         "procedure_datetime", "dx_mri_deltaT", "phecode", "description", "group",
                         "exclude_flag_AABTS", "age_at_scan_months")

#total number of rows across all the saved include tables
print(paste0("Total rows across all files in the list: ", sum(sapply(included_df_list, nrow))))

print(paste0("Binding rows based on columns specified and keeping only distinct rows."))
#merge all columns specified above and keep distinct rows
merged_df <- included_df_list %>%
  lapply(function(df){df[columns_to_merge]}) %>%
  bind_rows() %>% distinct()
print(paste0("Dimensions in the merged_df with distinct rows: ", dim(merged_df)))

merged_df$patID_phecode <- paste(merged_df$pat_id, merged_df$phecode, sep = "_")
print(paste0("Number of unique patID-phecode combinations: ", length(unique(merged_df$patID_phecode))))

#Get counts of all diagnoses under each patient
alldx_counts <- merged_df %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  #filter(phecode %in% phecodes_of_interest) %>%
  group_by(description) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

print(head(alldx_counts, 20))

#Get counts of neuropsych dx of interest under each patient
nsdx_counts <- merged_df %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(description) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

print(nsdx_counts)

#Get counts of neuropsych diagnoses under each patient BEFORE procedure date
nsdx_counts_before <- merged_df %>%
  filter(dx_mri_deltaT <= 0) %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(description) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

print(nsdx_counts_before)

#Get counts of neuropsych diagnoses under each patient AFTER procedure date
nsdx_counts_after <- merged_df %>%
  filter(dx_mri_deltaT > 0) %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(description) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

print(nsdx_counts_after)

print(paste0("generating plot now"))
plot_df <- merged_df %>%
  distinct(patID_phecode,.keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest & as.numeric(gestational_age_num) <= 42 & as.numeric(gestational_age_num) >= 24 )
#Plot the neuropsych diagnoses by GAB
nsdx_gab_plot <-
ggplot(plot_df, 
       aes(x = description, y = as.numeric(gestational_age_num), color = description, fill = description)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(width = 0.1, fill = "gray", alpha = 0.3, shape = 21) +
  labs(title = "Distribution of NeuroPsych Dx by Gestational Age",
       x = NULL,
       y = "Gestational Age (weeks)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure_filename <- "/mnt/arcus/lab/users/kafadare/GAB_cohort/figures/neuropsych_GAB_violin.png"
print(paste0("saving figure now as: ", figure_filename))
ggsave(figure_filename, plot = nsdx_gab_plot)

plot_obj_filename <- "/mnt/arcus/lab/users/kafadare/GAB_cohort/figures/neuropsych_GAB_violin.rdata"
print(paste0("saving plot obj now as: ", plot_obj_filename))
save(nsdx_gab_plot, file = plot_obj_filename)

df_filename <- "/mnt/arcus/lab/users/kafadare/GAB_cohort/included/merged_included.rds"
print(paste0("saving merged df now as: ", df_filename))
saveRDS(merged_df, df_filename)

