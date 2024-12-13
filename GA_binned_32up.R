# LOOK AT NEUROPSYCH / NEURODEV diagnoses
# Phecodes of potential interest noted in Notion
# Neuropsych Phecodes of interest
# 
# Source: PheWas Catalog phecodes ICD10
# 
# 315 - Developmental delays and disorders (includes mental retardation
#   315.1 Learning disorder
#   315.2 Speech and language disorders
#   315.3 Mental retardation
#                                           
# 313 - Pervasive developmental disorders (includes vague emotional disorders and even nonorganic enuresis)
#   313.1 ADHD
#   313.2 tics & stuttering
#   313.3 Autism
#   
# 300 - Anxiety, phobic, dissociative disorders
#   300.1 Anxiety
#     300.11  GAD
#     300.12 Agoraphobia, social phobia, panic disorders
#     300.12 Phobias
#   300.3 Obsessive compulsive disorders
#   300.4 Dysthymic disorder (why not under mood idk?)
#   300.8 Acute reaction to stress
#   300.9 PTSD
#   
# 296 - Mood
#   296.1 Bipolar
#   296.22 MDD
# 
# 295 - Schizophrenia and other psychotic disorders
#   295.1 Schizophrenia
#   295.2 Paranoid disorder
#   295.3 Psychosis

require(bigrquery)
require(tidyr)
require(dplyr)
require(ggplot2)

phecodes_of_interest <- c(315.1, 315.2, 313.1, 313.3, 300.1, 300.3, 296, 296.1, 296.22, 295, 295.1, 295.2, 295.3)

###code to phecode etc
#pathtoICDcode '/home/alexanderba/arcus-basics/R/ICDexclusion_code_March2024/'
icd10tophecode <- read.csv('/home/alexanderba/arcus-basics/R/ICDexclusion_code_March2024/Phecode_map_v1_2_icd10cm_beta.csv') #downloaded Oct12 2023  from https://phewascatalog.org/phecodes_icd10
phecode_defs <- read.csv('/home/alexanderba/arcus-basics/R/ICDexclusion_code_March2024/phecode_definitions1.2.csv')

excluded_phecodes <- read.csv('/home/alexanderba/arcus/shared/filter-by-dx/dx-filters/phecodes_with_exclusion_TS_and_AAB_19April2024.csv')
excluded_phecodes <- excluded_phecodes[which(excluded_phecodes$exclude_or_include_AAB_TS == 'exclude'),]

bq_auth()
proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]

## GA 28 to 32 weeks
my_query <- 'SELECT * FROM lab.GA_phecodes_28to32_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_28to32 <- bq_table_download(results, page_size=2000)

#create a column for exclusion based on AAB_TS phecode exclusion list
GA_phecodes_28to32$excluded <- ifelse(GA_phecodes_28to32$phecode %in% excluded_phecodes$phecode, "exclude", "include")
#create column for patient ID (+) phecode pairs
GA_phecodes_28to32$patID_phecode <- paste(GA_phecodes_28to32$pat_id, GA_phecodes_28to32$phecode, sep = "_")

#number of unique pairs
length(unique(GA_phecodes_28to32$patID_phecode)) #64,213
#number of unique patient IDs
length(unique(GA_phecodes_28to32$pat_id))#1583

#df including patients who do not have any lifetime recorded diagnosis in the exclusion list, agnostic to the time of scan vs diagnosis
no_exclude_28to32 <- GA_phecodes_28to32 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id)
#446 IDs without any lifetime diagnosis exclusions based on phecode

#**same dx repeated for pat_ids, here selecting only distinct pat_id + phecode pairs, but not necessarily the "first entry" of a specific phecode(dx) for a specific pat_id.

#Count table by patient, filtering out duplicate patient-phecode pairs as explained above
exclude_count_28to32 <- GA_phecodes_28to32 %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(pat_id) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))

#Count table for diagnoses that are excluded (not including diagnoses that aren't excluded)
exclude_count_byDx_28to32 <- GA_phecodes_28to32 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(phecode_str) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% filter(exclude_count != 0)  %>% arrange(desc(exclude_count))

#Looking at the distribution of how many excluded diagnoses per patient
density(exclude_count_28to32$exclude_count)

#PLOTS

#Plot histogram of excluded dx counts by patient (phecode-based)
hist(exclude_count_28to32$exclude_count,
     main = "28 to 32 weeks",
     xlab = "Number of Excludes",
     col = "blue",
     border = "black",
     breaks = 100)

#Plot the top 20 diagnoses that are excluded
ggplot(head(exclude_count_byDx_28to32, 20), aes(x = reorder(phecode_str, -exclude_count), y = exclude_count, fill = reorder(phecode_str, -exclude_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Excluded Diagnoses in 28 to 32 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# DIAGNOSES OF INTEREST

#Count table for diagnoses of interest defined at the top of this script, without any exclusions from the sample
neuropsychDx_count_28to32_noExclusions <- GA_phecodes_28to32 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

no_exclude_28to32 <- GA_phecodes_28to32 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id, .keep_all = TRUE)
dim(no_exclude_28to32)

#Count table for diagnoses of interest above, after lifetime exclusions from the sample
neuropsychDx_count_28to32_afterExclusions <- no_exclude_28to32 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

#Plot the diagnoses

#Plot the diagnoses of interest
ggplot(neuropsychDx_count_28to32_noExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis Before Exclusions, 28 to 32 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggplot(neuropsychDx_count_28to32_afterExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis After Exclusions, 28 to 32 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

rm(no_exclude_28to32)

rm(GA_phecodes_28to32)   


##32 to 34 weeks

my_query <- 'SELECT * FROM lab.GA_phecodes_32to34_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_32to34 <- bq_table_download(results, page_size=2000)

#create a column for exclusion based on AAB_TS phecode exclusion list
GA_phecodes_32to34$excluded <- ifelse(GA_phecodes_32to34$phecode %in% excluded_phecodes$phecode, "exclude", "include")
#create column for patient ID (+) phecode pairs
GA_phecodes_32to34$patID_phecode <- paste(GA_phecodes_32to34$pat_id, GA_phecodes_32to34$phecode, sep = "_")

#number of unique pairs
length(unique(GA_phecodes_32to34$patID_phecode)) #60,116
#number of unique patient IDs
length(unique(GA_phecodes_32to34$pat_id))#1566

#df including patients who do not have any lifetime recorded diagnosis in the exclusion list, agnostic to the time of scan vs diagnosis
no_exclude_32to34 <- GA_phecodes_32to34 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id)
dim(no_exclude_32to34)
#482 IDs without any lifetime diagnosis exclusions based on phecode

#**same dx repeated for pat_ids, here selecting only distinct pat_id + phecode pairs, but not necessarily the "first entry" of a specific phecode(dx) for a specific pat_id.

#Count table by patient, filtering out duplicate patient-phecode pairs as explained above
exclude_count_32to34 <- GA_phecodes_32to34 %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(pat_id) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))

#Count table for diagnoses that are excluded (not including diagnoses that aren't excluded)
exclude_count_byDx_32to34 <- GA_phecodes_32to34 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(phecode_str) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% filter(exclude_count != 0)  %>% arrange(desc(exclude_count))

#Looking at the distribution of how many excluded diagnoses per patient
density(exclude_count_32to34$exclude_count)

#PLOTS

#Plot histogram of excluded dx counts by patient (phecode-based)
hist(exclude_count_32to34$exclude_count,
     main = "32 to 34 weeks",
     xlab = "Number of Excludes",
     col = "blue",
     border = "black",
     breaks = 100)

#Plot the top 20 diagnoses that are excluded
ggplot(head(exclude_count_byDx_32to34, 20), aes(x = reorder(phecode_str, -exclude_count), y = exclude_count, fill = reorder(phecode_str, -exclude_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Excluded Diagnoses in 32 to 34 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# DIAGNOSES OF INTEREST

#Count table for diagnoses of interest defined at the top of this script, without any exclusions from the sample
neuropsychDx_count_32to34_noExclusions <- GA_phecodes_32to34 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

no_exclude_32to34 <- GA_phecodes_32to34 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id, .keep_all = TRUE)
dim(no_exclude_32to34)

#Count table for diagnoses of interest above, after lifetime exclusions from the sample
neuropsychDx_count_32to34_afterExclusions <- no_exclude_32to34 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

#Plot the diagnoses

#Plot the diagnoses of interest
ggplot(neuropsychDx_count_32to34_noExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis Before Exclusions, 32 to 34 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggplot(neuropsychDx_count_32to34_afterExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis After Exclusions, 32 to 34 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

rm(no_exclude_32to34)

rm(GA_phecodes_32to34)  


## 34 to 37 weeks
my_query <- 'SELECT * FROM lab.GA_phecodes_34to37_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_34to37 <- bq_table_download(results, page_size=2000)

#create a column for exclusion based on AAB_TS phecode exclusion list
GA_phecodes_34to37$excluded <- ifelse(GA_phecodes_34to37$phecode %in% excluded_phecodes$phecode, "exclude", "include")
#create column for patient ID (+) phecode pairs
GA_phecodes_34to37$patID_phecode <- paste(GA_phecodes_34to37$pat_id, GA_phecodes_34to37$phecode, sep = "_")

#number of unique pairs
length(unique(GA_phecodes_34to37$patID_phecode)) #213,792
#number of unique patient IDs
length(unique(GA_phecodes_34to37$pat_id))#5936

#df including patients who do not have any lifetime recorded diagnosis in the exclusion list, agnostic to the time of scan vs diagnosis
no_exclude_34to37 <- GA_phecodes_34to37 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id)
dim(no_exclude_34to37)
#1991 IDs without any lifetime diagnosis exclusions based on phecode

#**same dx repeated for pat_ids, here selecting only distinct pat_id + phecode pairs, but not necessarily the "first entry" of a specific phecode(dx) for a specific pat_id.

#Count table by patient, filtering out duplicate patient-phecode pairs as explained above
exclude_count_34to37 <- GA_phecodes_34to37 %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(pat_id) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))

#Count table for diagnoses that are excluded (not including diagnoses that aren't excluded)
exclude_count_byDx_34to37 <- GA_phecodes_34to37 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(phecode_str) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% filter(exclude_count != 0)  %>% arrange(desc(exclude_count))

#Looking at the distribution of how many excluded diagnoses per patient
density(exclude_count_34to37$exclude_count)

#PLOTS

#Plot histogram of excluded dx counts by patient (phecode-based)
hist(exclude_count_34to37$exclude_count,
     main = "34 to 37 weeks",
     xlab = "Number of Excludes",
     col = "blue",
     border = "black",
     breaks = 100)

#Plot the top 20 diagnoses that are excluded
ggplot(head(exclude_count_byDx_34to37, 20), aes(x = reorder(phecode_str, -exclude_count), y = exclude_count, fill = reorder(phecode_str, -exclude_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Excluded Diagnoses in 34 to 37 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# DIAGNOSES OF INTEREST

#Count table for diagnoses of interest defined at the top of this script, without any exclusions from the sample
neuropsychDx_count_34to37_noExclusions <- GA_phecodes_34to37 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

no_exclude_34to37 <- GA_phecodes_34to37 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id, .keep_all = TRUE)
dim(no_exclude_34to37)

#Count table for diagnoses of interest above, after lifetime exclusions from the sample
neuropsychDx_count_34to37_afterExclusions <- no_exclude_34to37 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

#Plot the diagnoses

#Plot the diagnoses of interest
ggplot(neuropsychDx_count_34to37_noExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis Before Exclusions, 34 to 37 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggplot(neuropsychDx_count_34to37_afterExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis After Exclusions, 34 to 37 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

rm(no_exclude_34to37)

rm(GA_phecodes_34to37)  


## 37 to 39 weeks
my_query <- 'SELECT * FROM lab.GA_phecodes_37to39_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_37to39 <- bq_table_download(results, page_size=2000)

#create a column for exclusion based on AAB_TS phecode exclusion list
GA_phecodes_37to39$excluded <- ifelse(GA_phecodes_37to39$phecode %in% excluded_phecodes$phecode, "exclude", "include")
#create column for patient ID (+) phecode pairs
GA_phecodes_37to39$patID_phecode <- paste(GA_phecodes_37to39$pat_id, GA_phecodes_37to39$phecode, sep = "_")

#number of unique pairs
length(unique(GA_phecodes_37to39$patID_phecode)) #350,380
#number of unique patient IDs
length(unique(GA_phecodes_37to39$pat_id))#10,133

#df including patients who do not have any lifetime recorded diagnosis in the exclusion list, agnostic to the time of scan vs diagnosis
no_exclude_37to39 <- GA_phecodes_37to39 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id)
dim(no_exclude_37to39)
#3478 IDs without any lifetime diagnosis exclusions based on phecode

#**same dx repeated for pat_ids, here selecting only distinct pat_id + phecode pairs, but not necessarily the "first entry" of a specific phecode(dx) for a specific pat_id.

#Count table by patient, filtering out duplicate patient-phecode pairs as explained above
exclude_count_37to39 <- GA_phecodes_37to39 %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(pat_id) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))

#Count table for diagnoses that are excluded (not including diagnoses that aren't excluded)
exclude_count_byDx_37to39 <- GA_phecodes_37to39 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(phecode_str) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% filter(exclude_count != 0)  %>% arrange(desc(exclude_count))

#Looking at the distribution of how many excluded diagnoses per patient
density(exclude_count_37to39$exclude_count)

#PLOTS

#Plot histogram of excluded dx counts by patient (phecode-based)
hist(exclude_count_37to39$exclude_count,
     main = "37 to 39 weeks",
     xlab = "Number of Excludes",
     col = "blue",
     border = "black",
     breaks = 100)

#Plot the top 20 diagnoses that are excluded
ggplot(head(exclude_count_byDx_37to39, 20), aes(x = reorder(phecode_str, -exclude_count), y = exclude_count, fill = reorder(phecode_str, -exclude_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Excluded Diagnoses in 37 to 39 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# DIAGNOSES OF INTEREST

#Count table for diagnoses of interest defined at the top of this script, without any exclusions from the sample
neuropsychDx_count_37to39_noExclusions <- GA_phecodes_37to39 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

no_exclude_37to39 <- GA_phecodes_37to39 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id, .keep_all = TRUE)
dim(no_exclude_37to39)

#Count table for diagnoses of interest above, after lifetime exclusions from the sample
neuropsychDx_count_37to39_afterExclusions <- no_exclude_37to39 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

#Plot the diagnoses

#Plot the diagnoses of interest
ggplot(neuropsychDx_count_37to39_noExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis Before Exclusions, 37 to 39 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggplot(neuropsychDx_count_37to39_afterExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis After Exclusions,37 to 39 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

rm(no_exclude_37to39)

rm(GA_phecodes_37to39)  

##39 to 41 weeks
my_query <- 'SELECT * FROM lab.GA_phecodes_39to41_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_39to41 <- bq_table_download(results, page_size=2000)

#create a column for exclusion based on AAB_TS phecode exclusion list
GA_phecodes_39to41$excluded <- ifelse(GA_phecodes_39to41$phecode %in% excluded_phecodes$phecode, "exclude", "include")
#create column for patient ID (+) phecode pairs
GA_phecodes_39to41$patID_phecode <- paste(GA_phecodes_39to41$pat_id, GA_phecodes_39to41$phecode, sep = "_")

#number of unique pairs
length(unique(GA_phecodes_39to41$patID_phecode)) #739,504
#number of unique patient IDs
length(unique(GA_phecodes_39to41$pat_id))#23737

#df including patients who do not have any lifetime recorded diagnosis in the exclusion list, agnostic to the time of scan vs diagnosis
no_exclude_39to41 <- GA_phecodes_39to41 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id)
dim(no_exclude_39to41)
#9792 IDs without any lifetime diagnosis exclusions based on phecode

#**same dx repeated for pat_ids, here selecting only distinct pat_id + phecode pairs, but not necessarily the "first entry" of a specific phecode(dx) for a specific pat_id.

#Count table by patient, filtering out duplicate patient-phecode pairs as explained above
exclude_count_39to41 <- GA_phecodes_39to41 %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(pat_id) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))

#Count table for diagnoses that are excluded (not including diagnoses that aren't excluded)
exclude_count_byDx_39to41 <- GA_phecodes_39to41 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(phecode_str) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% filter(exclude_count != 0)  %>% arrange(desc(exclude_count))

#Looking at the distribution of how many excluded diagnoses per patient
density(exclude_count_39to41$exclude_count)

#PLOTS

#Plot histogram of excluded dx counts by patient (phecode-based)
hist(exclude_count_39to41$exclude_count,
     main = "39 to 41 weeks",
     xlab = "Number of Excludes",
     col = "blue",
     border = "black",
     breaks = 100)

#Plot the top 20 diagnoses that are excluded
ggplot(head(exclude_count_byDx_39to41, 20), aes(x = reorder(phecode_str, -exclude_count), y = exclude_count, fill = reorder(phecode_str, -exclude_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Excluded Diagnoses in 39 to 41 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# DIAGNOSES OF INTEREST

#Count table for diagnoses of interest defined at the top of this script, without any exclusions from the sample
neuropsychDx_count_39to41_noExclusions <- GA_phecodes_39to41 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

no_exclude_39to41 <- GA_phecodes_39to41 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id, .keep_all = TRUE)
dim(no_exclude_39to41)

#Count table for diagnoses of interest above, after lifetime exclusions from the sample
neuropsychDx_count_39to41_afterExclusions <- no_exclude_39to41 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

#Plot the diagnoses

#Plot the diagnoses of interest
ggplot(neuropsychDx_count_39to41_noExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis Before Exclusions, 39 to 41 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggplot(neuropsychDx_count_39to41_afterExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis After Exclusions,39 to 41 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Comorbidity of diagnoses of interest with the "exclusion" diagnoses

# Create a binary presence table for each patient with the relevant diagnosis (dx of interest and dx of exclusion)
diagnosis_table <- GA_phecodes_39to41 %>%
  filter(phecode %in% c(phecodes_of_interest, excluded_phecodes$phecode)) %>%
  select(pat_id, phecode_str) %>%
  pivot_wider(names_from = phecode_str, values_from = phecode_str, values_fill = 0, values_fn = length) %>%
  mutate(across(-pat_id, ~ ifelse(. > 0, 1, 0)))

# Calculate correlations
cor_matrix <- diagnosis_table %>%
  select(-pat_id) %>%
  cor()

# Extract and order pairs
phecodes_of_interest_names <- icd10tophecode[which(icd10tophecode$phecode %in% phecodes_of_interest), 'phecode_str']

cor_pairs <- as.data.frame(as.table(cor_matrix)) %>%
  filter(Var1 %in% phecodes_of_interest_names | Var2 %in% phecodes_of_interest_names) %>%
  filter(!(Var1 %in% phecodes_of_interest_names & Var2 %in% phecodes_of_interest_names)) %>%
  #filter(!(Var1 == Var2)) %>%
  arrange(desc(Freq))

# Calculate co-occurrence counts: within summarise fct, list applies the custom function to each column respectively, here is does an element-wise multiplication to count the number of co-occurences between diagnoses
## !!! The ph_excl column is just a repeat of the ph_interest column !!!

cooccurrence_counts <- diagnosis_table %>%
  summarise(across(names(diagnosis_table)[which(names(diagnosis_table) %in% phecodes_of_interest_names)], 
                   list(~ colSums(. * select(diagnosis_table, names(diagnosis_table)[which(names(diagnosis_table) %in% excluded_phecodes$phenotype)]))),
                   .names = "cooccur_{.col}_{col}")) %>%
  pivot_longer(cols = starts_with("cooccur"), 
               names_to = c("ph_interest", "ph_excl"), 
               names_pattern = "cooccur_(.*)_(.*)", 
               values_to = "count") %>%
  arrange(desc(count))

cooccurrence_counts

head(cor_pairs, 40)

rm(no_exclude_39to41)

rm(GA_phecodes_39to41)

##Over 40 weeks

my_query <- 'SELECT * FROM lab.GA_phecodes_over40_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_over40 <- bq_table_download(results, page_size=2000)

#create a column for exclusion based on AAB_TS phecode exclusion list
GA_phecodes_over40$excluded <- ifelse(GA_phecodes_over40$phecode %in% excluded_phecodes$phecode, "exclude", "include")
#create column for patient ID (+) phecode pairs
GA_phecodes_over40$patID_phecode <- paste(GA_phecodes_over40$pat_id, GA_phecodes_over40$phecode, sep = "_")
#number of unique pairs
length(unique(GA_phecodes_over40$patID_phecode)) #100,326
#number of unique patient IDs
length(unique(GA_phecodes_over40$pat_id))#3186

#df including patients who do not have any lifetime recorded diagnosis in the exclusion list, agnostic to the time of scan vs diagnosis
no_exclude_over40 <- GA_phecodes_over40 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id)
dim(no_exclude_over40)
#1302 IDs without any lifetime diagnosis exclusions based on phecode

#**same dx repeated for pat_ids, here selecting only distinct pat_id + phecode pairs, but not necessarily the "first entry" of a specific phecode(dx) for a specific pat_id.

#Count table by patient, filtering out duplicate patient-phecode pairs as explained above
exclude_count_over40 <- GA_phecodes_over40 %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(pat_id) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))

#Count table for diagnoses that are excluded (not including diagnoses that aren't excluded)
exclude_count_byDx_over40 <- GA_phecodes_over40 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(phecode_str) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% filter(exclude_count != 0)  %>% arrange(desc(exclude_count))

#Looking at the distribution of how many excluded diagnoses per patient
density(exclude_count_over40$exclude_count)

#PLOTS

#Plot histogram of excluded dx counts by patient (phecode-based)
hist(exclude_count_over40$exclude_count,
     main = "41 weeks and over",
     xlab = "Number of Excludes",
     col = "blue",
     border = "black",
     breaks = 100)

#Plot the top 20 diagnoses that are excluded
ggplot(head(exclude_count_byDx_over40, 20), aes(x = reorder(phecode_str, -exclude_count), y = exclude_count, fill = reorder(phecode_str, -exclude_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Excluded Diagnoses in 41 weeks and over",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# DIAGNOSES OF INTEREST

#Count table for diagnoses of interest defined at the top of this script, without any exclusions from the sample
neuropsychDx_count_over40_noExclusions <- GA_phecodes_over40 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

no_exclude_over40 <- GA_phecodes_over40 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id, .keep_all = TRUE)
dim(no_exclude_over40)

#Count table for diagnoses of interest above, after lifetime exclusions from the sample
neuropsychDx_count_over40_afterExclusions <- no_exclude_over40 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

#Plot the diagnoses

#Plot the diagnoses of interest
ggplot(neuropsychDx_count_over40_noExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis Before Exclusions, 40 weeks and over",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggplot(neuropsychDx_count_over40_afterExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis After Exclusions, 40 weeks and over",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

rm(no_exclude_over40)

rm(GA_phecodes_over40)  


