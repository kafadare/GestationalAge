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

bq_auth()
proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]

query_inc <- glue("select * from lab.GA_procedure_orders_EK where proc_ord_id NOT IN (select cast(proc_ord_id as string) from lab.proc_ord_id_exclude_masterlist_EK) AND proc_ord_id IN (select proc_ord_id from lab.proc_ids_inDxTables_EK);")

print(paste0("submitting query: ", query_inc))

df_inc <- bq_project_query(proj_id, query_inc) %>% bq_table_download(., page_size=3000)

print(paste0("Included procedures downloaded, dimensions: ", dim(df_inc)))

print(paste0("Number of unique proc_ord_ids: ", length(unique(df_inc$proc_ord_id))))

print(paste0("Number of unique pat_ids: ", length(unique(df_inc$pat_id))))

query_all <- glue("select * from lab.GA_procedure_orders_EK where proc_ord_id IN (select proc_ord_id from lab.proc_ids_inDxTables_EK);")

print(paste0("submitting query: ", query_all))

df_all <- bq_project_query(proj_id, query_all) %>% bq_table_download(., page_size=3000)

print(paste0("All procedures downloaded, dimensions: ", dim(df_all)))

print(paste0("Number of unique proc_ord_ids: ", length(unique(df_all$proc_ord_id))))

print(paste0("Number of unique pat_ids: ", length(unique(df_all$pat_id))))

#define procedure start time
df_inc$procedure_datetime <- df_inc$start_datetime
df_inc[is.na(df_inc$procedure_datetime),'procedure_datetime'] <- df_inc[is.na(df_inc$procedure_datetime), 'proc_ord_datetime']
print(paste0("How many procedure_datetime is NA: ",sum(is.na(df_inc$procedure_datetime))))#0

df_inc$age_at_scan_months <- as.numeric(df_inc$procedure_datetime - df_inc$dob)/(60*60*24*30.43)#in months

df_all$procedure_datetime <- df_all$start_datetime
df_all[is.na(df_all$procedure_datetime),'procedure_datetime'] <- df_all[is.na(df_all$procedure_datetime), 'proc_ord_datetime']
print(paste0("How many procedure_datetime is NA: ",sum(is.na(df_all$procedure_datetime))))#0

df_all$age_at_scan_months <- as.numeric(df_all$procedure_datetime - df_all$dob)/(60*60*24*30.43)#in months


#Look at distribution of GAB/age before & after exclusion of exclusion phecodes by date dx < date MRI, buffer 30 days.
  GAB_Age_plot <- 
    df_all %>% 
    distinct(proc_ord_id, .keep_all = TRUE) %>%
    ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
    labs(
      title = "Pre-Exclusion Sample, N = 95,315 (proc IDs)",
      x = "Age At Scan (months)",
      y = "GAB (weeks)" 
    ) +
    scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))
  
  #look at distribution plot again
  GAB_Age_excluded_plot <- 
    df_inc %>%
    distinct(proc_ord_id, .keep_all = TRUE) %>%
    ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
    labs(
      title = "Post-Exclusion Sample, N = 28690 (proc IDs)",
      x = "Age At Scan (months)",
      y = "GAB (weeks)" 
    ) +
    scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))

  
#Look at the grading distribution of the included proc_ord_ids in this cohort
query_grades <- glue("SELECT proc_ord_id, grader_name, grade, grade_category, report_origin_table, grade_date FROM lab.grader_table_with_metadata_project_independent;")
  
print(paste0("submitting query: ", query_grades))
  
grades <- bq_project_query(proj_id, query_grades) %>% bq_table_download(., page_size=3000)
  
print(paste0("Grader table, dimensions: ", dim(grades)))

print(paste0("Number of unique proc_ord_ids: ", length(unique(grades$proc_ord_id))))

head(grades,20)
table(grades$grader_name)
table(grades$grade_category)

#Merge included proc_ord_ids with grader table
df_inc_graded <- merge(df_inc, grades, by = "proc_ord_id")
print(paste0("Merged, dimensions: ", dim(df_inc_graded)))
print(paste0("Merged, number of unique proc_ord_ids: ", length(unique(df_inc_graded$proc_ord_id))))
print(paste0("Merged, number of unique patient ids: ", length(unique(df_inc_graded$pat_id))))

df_inc_graded <- df_inc_graded %>% 
  subset(grade >= 0 & grade <= 2) %>%
  group_by(proc_ord_id) %>%
  mutate(mean_grade = mean(grade))

range(df_inc_graded$mean_grade) #0 to 2

print(paste0("Remove absurd grades, number of unique proc_ord_ids: ", length(unique(df_inc_graded$proc_ord_id))))

GAB_grades_histogram <- 
  df_inc_graded %>% distinct(proc_ord_id, .keep_all = TRUE) %>%
  ggplot(., aes(x = gestational_age_num, fill = mean_grade)) + 
  geom_histogram(binwidth = 1, position = "stack", stat = "count") +
  labs(
    title = "Distribution of Mean Grade Across Gestational Age",
    x = "Age",
    y = "Count" 
  ) +
  scale_fill_viridis_c(option = "viridis") +
  theme_minimal()

GAB_Age_grades_heatmap <- 
  df_inc_graded %>% distinct(proc_ord_id, .keep_all = TRUE) %>%
  ggplot(., aes(x = age_at_scan_months, y = gestational_age_num, color = mean_grade, fill = mean_grade)) +
  geom_tile() +
  #geom_bin2d(bins = 100)+
  labs(
    title = "HeatMap of Graded Reports by Age and GAB",
    x = "Age At Scan (months)",
    y = "GAB (weeks)" 
  ) +
  scale_fill_viridis_c(option = "viridis") +
  scale_color_viridis_c(option = "viridis") +
  theme_minimal()

dim(df_inc_graded %>% filter(mean_grade >= 1.5) %>% distinct(proc_ord_id))
dim(df_inc_graded %>% filter(mean_grade >= 1.5) %>% distinct(pat_id))

GAB_Age_normalGrade_heatmap <- 
df_inc_graded %>%
  distinct(proc_ord_id, .keep_all = TRUE) %>%
  filter(mean_grade >= 1.5) %>%
  ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
  labs(
    title = "Mean Grade >= 1.5, N = 3258 (proc IDs)",
    x = "Age At Scan (months)",
    y = "GAB (weeks)" 
  ) +
  scale_fill_viridis_c(option = "viridis")

dim(df_inc_graded %>% filter(mean_grade <= 0.5) %>% distinct(proc_ord_id))
dim(df_inc_graded %>% filter(mean_grade <= 0.5) %>% distinct(pat_id))

GAB_Age_abnormalGrade_heatmap <- 
  df_inc_graded %>%
  distinct(proc_ord_id, .keep_all = TRUE) %>%
  filter(mean_grade <= 0.5) %>%
  ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
  labs(
    title = "Mean Grade <= 0.5, N = 3196 (proc IDs)",
    x = "Age At Scan (months)",
    y = "GAB (weeks)" 
  ) +
  scale_fill_viridis_c(option = "viridis")

