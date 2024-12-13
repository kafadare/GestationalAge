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

bq_auth()
proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]
rowcount_query <- 'SELECT count(*) FROM lab.proc_encDx_filter_EK'
rowcount <- bq_project_query(proj_id, rowcount_query) %>% bq_table_download()

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
}

print(paste0("NUMBER OF CHUNKS is: ", chunks))

size <- ceiling(rowcount/chunks)
print(paste0("SIZE is: ", size))

saveRDS <- FALSE

plot <- FALSE

upload <- FALSE

offset <- n-1

query <- glue("SELECT * FROM lab.proc_encDx_filter_EK LIMIT {size} OFFSET {size * offset}")

print(paste0("submitting query: ", query))

encDx <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500)

print(paste0("encDx downloaded, dimensions: ", dim(encDx)))

#define procedure start time
encDx$procedure_datetime <- encDx$start_datetime
encDx[is.na(encDx$procedure_datetime),'procedure_datetime'] <- encDx[is.na(encDx$procedure_datetime), 'proc_ord_datetime']
print(paste0("How many procedure_datetime is NA: ",sum(is.na(encDx$procedure_datetime))))#0

#define diagnosis datetime
encDx$diagnosis_datetime <- encDx$effective_datetime

#calculate dx_mri_deltaT
encDx$dx_mri_deltaT <- (encDx$procedure_datetime - encDx$diagnosis_datetime)/(60*60*24)
print(paste0("How many dx_mri_deltaT is NA: ",sum(is.na(encDx$dx_mri_deltaT))))#0

#hist(as.numeric(encDx$dx_mri_deltaT), breaks = 100)

#Map onto Phecodes

library(PheWAS)
print(paste0("Starting the icd_code - PheWAS conversion"))

encDx$ind <- 1:nrow(encDx)#create a column for row index to match the phecodes later.

##Splitting the ICD code and unnesting -- there are multiple ICD codes per cell divided by commas.

encDx_icd10 <- encDx %>%
  add_column(vocabulary_id = "ICD10CM",
             code = encDx$icd10_list) %>%
  select(ind, vocabulary_id, code) %>%
  drop_na() %>%
  distinct() %>%
  mutate(code = str_split(code, ', ')) %>%
  unnest_longer(code)%>%
  mapCodesToPhecodes() %>%
  distinct(.keep_all = T)

encDx_icd9 <- encDx %>%
  add_column(vocabulary_id = "ICD9CM",
             code = encDx$icd9_list) %>%
  select(ind, vocabulary_id, code) %>%
  drop_na() %>%
  distinct() %>%
  mutate(code = str_split(code, ', ')) %>%
  unnest_longer(code)%>%
  mapCodesToPhecodes() %>%
  distinct(.keep_all = T)

encDx_icd <- encDx_icd10 %>%
  bind_rows(encDx_icd9) %>%
  distinct(.keep_all = T)

rm(encDx_icd9)
rm(encDx_icd10)

# Merge phecodes with original dx, remove duplicates
print(paste0("Merging encDx with the converted phecodes"))

encDx <- encDx %>%
  merge(encDx_icd, by = "ind", all = T) %>%
  addPhecodeInfo()

encDx$exclude_flag_AABTS <- ifelse(encDx$phecode %in% excluded_phecodes$phecode, TRUE, FALSE)

if(saveRDS == TRUE){saveRDS(encDx, "/mnt/arcus/lab/users/kafadare/encDx_part1_phecodes.rds")}

encDx$age_at_scan_months <- as.numeric(encDx$procedure_datetime - encDx$dob)/(60*60*24*30.43)#in months

#exclude patient_scan if they have an exclusion phecode diagnosis within 30 days or before MRI procedure date.
print(paste0("Number of unique procedure IDs in this split: ", length(unique(encDx$proc_ord_id))))#
proc_ord_ids_excluded <- encDx %>% 
  filter((exclude_flag_AABTS == TRUE & dx_mri_deltaT > -30)) %>%
  distinct(proc_ord_id)
print(paste0("Number of unique procedure IDs that are excluded in this split: ", length(unique(proc_ord_ids_excluded))))#

#two column df for all distinct procedure order ids with an include/exclude flag
proc_ord_ids_all <- distinct(encDx,proc_ord_id)
proc_ord_ids_all$exclude <- ifelse(proc_ord_ids_all$proc_ord_id %in% proc_ord_ids_excluded$proc_ord_id, TRUE, FALSE)
print(paste0("Number of unique procedure IDs that have EXclude flag: ", sum(proc_ord_ids_all$exclude == TRUE)))#
print(paste0("Number of unique procedure IDs that have INclude flag: ", sum(proc_ord_ids_all$exclude == FALSE)))#
print(paste0("Number of unique procedure IDs that have have NA in exclude flas: ", sum(is.na(proc_ord_ids_all$exclude))))#

procedures_excluded <- encDx %>% 
  subset(!(proc_ord_id %in% proc_ord_ids_excluded))
print(paste0("Dimensions of procedures_excluded df: ", dim(procedures_excluded)))#

procedures_included <- encDx %>%
  subset(!(proc_ord_id %in% proc_ord_ids_excluded))
print(paste0("Dimensions of procedures_included df: ", dim(procedures_included)))#

filename_ids <- glue("/mnt/arcus/lab/users/kafadare/GAB_cohort/excluded_procIDs/proc_ord_ids_encDx_part{n}.csv")
print(paste0("Writing file: ", filename_ids))
saveRDS(proc_ord_ids_all, filename_ids)

filename_ex <- glue("/mnt/arcus/lab/users/kafadare/GAB_cohort/excluded/procedures_exclude_encDx_part{n}.rds")
print(paste0("Writing file: ", filename_ex))
saveRDS(procedures_excluded, filename_ex)

filename_inc <- glue("/mnt/arcus/lab/users/kafadare/GAB_cohort/included/procedures_include_encDx_part{n}.rds")
print(paste0("Writing file: ", filename_inc))
saveRDS(procedures_included, filename_inc)

if (upload == TRUE){
#Upload phecodes_excluded to SQL arcus
#procedures_excluded <- read.csv(filename)
bq_auth()
#proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]
print(paste0("writing into SQL table now."))

bq_con <- dbConnect(bigrquery::bigquery(), project = proj_id)
dest_table = bq_table(project = proj_id, dataset = "lab", table = "proc_ord_id_exclude_encDx_EK")
bq_table_upload(
  x=dest_table, 
  values= proc_ord_ids_excluded, 
  create_disposition='CREATE_IF_NEEDED', 
  write_disposition='WRITE_APPEND')
}
#Look at distribution of GAB/age before & after exclusion of exclusion phecodes by date dx < date MRI, buffer 30 days.
if(plot == TRUE){
GAB_Age_plot <- 
  encDx %>% distinct(proc_ord_id, .keep_all = TRUE) %>%
  ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
  scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))

#look at distribution plot again
GAB_Age_excluded_plot <- 
  encDx %>% 
  filter(!(proc_ord_id %in% procedures_excluded)) %>%
  distinct(proc_ord_id, .keep_all = TRUE) %>%
  ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
  scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))
}
