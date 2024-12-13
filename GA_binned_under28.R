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

####Under 28 weeks

my_query <- 'SELECT * FROM lab.GA_phecodes_under28_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_under28 <- bq_table_download(results, page_size=2000)

#create a column for exclusion based on AAB_TS phecode exclusion list
GA_phecodes_under28$excluded <- ifelse(GA_phecodes_under28$phecode %in% excluded_phecodes$phecode, "exclude", "include")
#create column for patient ID (+) phecode pairs
GA_phecodes_under28$patID_phecode <- paste(GA_phecodes_under28$pat_id, GA_phecodes_under28$phecode, sep = "_")

#number of unique pairs
length(unique(GA_phecodes_under28$patID_phecode)) #72,517
#number of unique patient IDs
length(unique(GA_phecodes_under28$pat_id))#1404

#df including patients who do not have any lifetime recorded diagnosis in the exclusion list, agnostic to the time of scan vs diagnosis
no_exclude_under28 <- GA_phecodes_under28 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id)
#247 IDs without any lifetime diagnosis exclusions based on phecode

####### This chunk is just trying out counts absolute (without filtering for distinct patient_id/phecode pairs), better to do it with the filtering going forwards as more informative. Can do it this way to also include the information about the # of visits/admissions.

# exclude_count_under28 <- GA_phecodes_under28 %>%
#   group_by(pat_id) %>%
#   summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))
# 
# #exploring the number of diagnoses and exclusions ...
# range(exclude_count_under28$exclude_count) # (0,5362)
# exclude_count_under28[which(exclude_count_under28$exclude_count == 5362),]
# # !!! HMRWUYN6 with 5362 "exclusions"!!!
# View(GA_phecodes_under28[which(GA_phecodes_under28$pat_id == "HMRWUYN6"),]) #Turns out nothing super weird, just a patient with a bunch of diagnoses/visits. There are some patients like this on the list, seems like most are chronic disease patients with lots of visits-admissions/diagnosis entries - which makes sense.

####### Chunk over

#**same dx repeated for pat_ids, here selecting only distinct pat_id + phecode pairs, but not necessarily the "first entry" of a specific phecode(dx) for a specific pat_id.

#Count table by patient, filtering out duplicate patient-phecode pairs as explained above
exclude_count_under28 <- GA_phecodes_under28 %>%
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(pat_id) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% arrange(desc(exclude_count))

#Count table for diagnoses that are excluded (not including diagnoses that aren't excluded)
exclude_count_byDx_under28 <- GA_phecodes_under28 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  group_by(phecode_str) %>%
  summarize(exclude_count = sum(excluded == "exclude")) %>% filter(exclude_count != 0)  %>% arrange(desc(exclude_count))

#Looking at the distribution of how many excluded diagnoses per patient
density(exclude_count_under28$exclude_count)
hist(exclude_count_under28$exclude_count,
     main = "Under 28 weeks",
     xlab = "Number of Excludes",
     col = "blue",
     border = "black",
     breaks = 100)

#Plot the top 20 diagnoses that are excluded
ggplot(head(exclude_count_byDx_under28, 20), aes(x = reorder(phecode_str, -exclude_count), y = exclude_count, fill = reorder(phecode_str, -exclude_count))) +
         geom_bar(stat = "identity", position = position_dodge()) +
         labs(title = "Incidence of Excluded Diagnoses in Under 28 Weeks",
              x = "Diagnosis",
              y = "Incidence in Sample") +
         theme_minimal() +
         theme(axis.text.x = element_blank())


# DIAGNOSES OF INTEREST

#Count table for diagnoses of interest defined at the top of this script, without any exclusions from the sample
neuropsychDx_count_under28_noExclusions <- GA_phecodes_under28 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

no_exclude_under28 <- GA_phecodes_under28 %>%
  group_by(pat_id) %>%
  filter(!any(excluded == "exclude")) %>%
  distinct(pat_id, .keep_all = TRUE)
dim(no_exclude_under28)

#Count table for diagnoses of interest above, after lifetime exclusions from the sample
neuropsychDx_count_under28_afterExclusions <- no_exclude_under28 %>% 
  distinct(patID_phecode, .keep_all = TRUE) %>%
  filter(phecode %in% phecodes_of_interest) %>%
  group_by(phecode_str) %>%
  summarize(diagnosis_count = n()) %>% 
  arrange(desc(diagnosis_count))

#Plot the diagnoses

#Plot the diagnoses of interest
ggplot(neuropsychDx_count_under28_noExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis Before Exclusions, Under 28 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggplot(neuropsychDx_count_under28_afterExclusions, aes(x = reorder(phecode_str, -diagnosis_count), y = diagnosis_count, fill = reorder(phecode_str, -diagnosis_count))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Incidence of Neuropsych Diagnosis After Exclusions, Under 28 weeks",
       x = "Diagnosis",
       y = "Incidence in Sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

rm(no_exclude_under28)
rm(GA_phecodes_under28)  

## PHECODES & PROCEDURE IDs --- EXCLUSION BY DATE
my_query <- 'SELECT * FROM lab.GA_phecodes_procedures_under28_EK'
results <- bq_project_query(proj_id, my_query)
GA_phecodes_procedures_under28 <- bq_table_download(results, page_size=2000)

#433098 unique
pat_phecode_encounter <- paste(GA_phecodes_procedures_under28$pat_id, 
                               GA_phecodes_procedures_under28$phecode, 
                               GA_phecodes_procedures_under28$encounter_id, sep = "_") %>% 
                        unique()

#68153 unique
pat_phecode <- paste(GA_phecodes_procedures_under28$pat_id, 
                               GA_phecodes_procedures_under28$phecode, sep = "_") %>% 
                        unique()

#encounters_dates <- GA_phecodes_procedures_under28 %>%
  #select(encounter_id, proc_ord_id, proc_ord_datetime, start_datetime, contact_date, pat_enc_date_real, effective_datetime, proc_ord_datetime, prob_entry_datetime, prob_noted_datetime)

#Procedure/Scan Date
sum(is.na(GA_phecodes_procedures_under28$proc_ord_datetime))
#0
sum(is.na(GA_phecodes_procedures_under28$start_datetime))
#334760
proc_start_datediff <- difftime(GA_phecodes_procedures_under28$start_datetime, GA_phecodes_procedures_under28$proc_ord_datetime, unit = "days")
#max difference in days is 8.64 days
#use start_datetime when available, proc_ord_datetime if not.

#Encounter/Dx Date -- Date associated with diagnosis
sum(is.na(GA_phecodes_procedures_under28$prob_noted_datetime))
#846091
sum(is.na(GA_phecodes_procedures_under28$prob_entry_datetime))
#839556
sum(is.na(GA_phecodes_procedures_under28$contact_date))
#0
sum(is.na(GA_phecodes_procedures_under28$pat_enc_date_real))
#0
sum(is.na(GA_phecodes_procedures_under28$effective_datetime))
#114

proc_noted_entry_datediff <- difftime(GA_phecodes_procedures_under28$prob_noted_datetime, GA_phecodes_procedures_under28$prob_entry_datetime, unit = "days")
range(proc_noted_entry_datediff, na.rm = T) #-7525, 0

proc_entry_contact_datediff <- difftime(GA_phecodes_procedures_under28$prob_entry_datetime, GA_phecodes_procedures_under28$contact_date, unit = "days")
range(proc_entry_contact_datediff, na.rm = T) #-77, 1612

proc_noted_contact_datediff <- difftime(GA_phecodes_procedures_under28$prob_noted_datetime, GA_phecodes_procedures_under28$contact_date, unit = "days")
range(proc_noted_contact_datediff, na.rm = T) #-7523, 1019

proc_noted_effective_datediff <- difftime(GA_phecodes_procedures_under28$prob_noted_datetime, GA_phecodes_procedures_under28$effective_datetime, unit = "days")
range(proc_noted_effective_datediff, na.rm = T) #-7523, 1017

cor.test(as.numeric(proc_noted_contact_datediff), as.numeric(proc_noted_effective_datediff)) # cor 0.9998
cor.test(as.numeric(proc_noted_contact_datediff), as.numeric(proc_noted_entry_datediff)) # cor 0.995

#Use prob_noted_datetime when available, otherwise use contact_date

GA_phecodes_procedures_under28 <- GA_phecodes_procedures_under28 %>%
  select(-c(effective_datetime, prob_entry_datetime, pat_enc_date_real)) %>% 
  mutate(age_at_scan_days = difftime(start_datetime, dob, unit = "days"), 
         age_at_dx_days = difftime(prob_noted_datetime, dob, unit = "days"))

sum(is.na(GA_phecodes_procedures_under28$age_at_scan_days)) #334760
sum(is.na(GA_phecodes_procedures_under28$age_at_dx_days)) #846091

test <- GA_phecodes_procedures_under28 %>%
  filter(!is.na(age_at_scan_days)) %>% 
  distinct(proc_ord_id)#1361/2011 proc_ids have age at scan when using start_datetime.

#fill out empty age at scan using proc ord date
GA_phecodes_procedures_under28[is.na(GA_phecodes_procedures_under28$age_at_scan_days),"age_at_scan_days"] <- difftime(GA_phecodes_procedures_under28[is.na(GA_phecodes_procedures_under28$age_at_scan_days),] $proc_ord_datetime, GA_phecodes_procedures_under28[is.na(GA_phecodes_procedures_under28$age_at_scan_days),]$dob, units = "days")
sum(is.na(GA_phecodes_procedures_under28$age_at_scan_days)) #0

#fill out empty age at dx using contact date
GA_phecodes_procedures_under28[is.na(GA_phecodes_procedures_under28$age_at_dx_days),"age_at_dx_days"] <- difftime(GA_phecodes_procedures_under28[is.na(GA_phecodes_procedures_under28$age_at_dx_days),] $contact_date, GA_phecodes_procedures_under28[is.na(GA_phecodes_procedures_under28$age_at_dx_days),]$dob, units = "days")
sum(is.na(GA_phecodes_procedures_under28$age_at_dx_days)) #0


##Want to see age distribution of scans, but filtering by distinct proc_id
hist(as.numeric(GA_phecodes_procedures_under28$age_at_scan_days))
