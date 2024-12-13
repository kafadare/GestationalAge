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
rowcount_query <- 'SELECT count(*) FROM lab.proc_procDx_filter_EK'
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

upload <- TRUE

offset <- n-1

query <- glue("SELECT * FROM lab.proc_procDx_filter_EK LIMIT {size} OFFSET {size * offset}")

print(paste0("submitting query: ", query))

procDx <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=2000)

print(paste0("procDx downloaded, dimensions: ", dim(procDx)))

#define procedure start time
procDx$procedure_datetime <- procDx$start_datetime
procDx[is.na(procDx$procedure_datetime),'procedure_datetime'] <- procDx[is.na(procDx$procedure_datetime), 'proc_ord_datetime']
print(paste0("How many procedure_datetime is NA: ",sum(is.na(procDx$procedure_datetime))))#0

#define diagnosis datetime
procDx$diagnosis_datetime <- procDx$effective_datetime

#calculate dx_mri_deltaT
procDx$dx_mri_deltaT <- (procDx$procedure_datetime - procDx$diagnosis_datetime)/(60*60*24)
print(paste0("How many dx_mri_deltaT is NA: ",sum(is.na(procDx$dx_mri_deltaT))))#0

#hist(as.numeric(procDx$dx_mri_deltaT), breaks = 100)

#Map onto Phecodes

library(PheWAS)
print(paste0("Starting the icd_code - PheWAS conversion"))

procDx$ind <- 1:nrow(procDx)#create a column for row index to match the phecodes later.

##Splitting the ICD code and unnesting -- there are multiple ICD codes per cell divided by commas.

procDx_icd10 <- procDx %>%
  add_column(vocabulary_id = "ICD10CM",
             code = procDx$icd10_list) %>%
  select(ind, vocabulary_id, code) %>%
  drop_na() %>%
  distinct() %>%
  mutate(code = str_split(code, ', ')) %>%
  unnest_longer(code)%>%
  mapCodesToPhecodes() %>%
  distinct(.keep_all = T)

procDx_icd9 <- procDx %>%
  add_column(vocabulary_id = "ICD9CM",
             code = procDx$icd9_list) %>%
  select(ind, vocabulary_id, code) %>%
  drop_na() %>%
  distinct() %>%
  mutate(code = str_split(code, ', ')) %>%
  unnest_longer(code)%>%
  mapCodesToPhecodes() %>%
  distinct(.keep_all = T)

procDx_icd <- procDx_icd10 %>%
  bind_rows(procDx_icd9) %>%
  distinct(.keep_all = T)

rm(procDx_icd9)
rm(procDx_icd10)

# Merge phecodes with original dx, remove duplicates
print(paste0("Merging procDx with the converted phecodes"))

procDx <- procDx %>%
  merge(procDx_icd, by = "ind", all = T) %>%
  addPhecodeInfo()

procDx$exclude_flag_AABTS <- ifelse(procDx$phecode %in% excluded_phecodes$phecode, TRUE, FALSE)

if(saveRDS == TRUE){saveRDS(procDx, "/mnt/arcus/lab/users/kafadare/procDx_part1_phecodes.rds")}

procDx$age_at_scan_months <- as.numeric(procDx$procedure_datetime - procDx$dob)/(60*60*24*30.43)#in months

#exclude patient_scan if they have an exclusion phecode diagnosis within 30 days or before MRI procedure date.
print(paste0("Number of unique procedure IDs in this split: ", length(unique(procDx$proc_ord_id))))#

proc_ord_ids_excluded <- procDx %>% 
  filter((exclude_flag_AABTS == TRUE & dx_mri_deltaT > -30)) %>%
  distinct(proc_ord_id)
print(paste0("Number of unique procedure IDs that are excluded in this split: ", length(unique(proc_ord_ids_excluded))))#

procedures_excluded <- procDx %>% 
  filter((exclude_flag_AABTS == TRUE & dx_mri_deltaT > -30))
print(paste0("Dimensions of procedures_excluded df: ", dim(procedures_excluded)))#

procedures_included <- procDx %>%
  subset(!(proc_ord_id %in% proc_ord_ids_excluded))
print(paste0("Dimensions of procedures_included df: ", dim(procedures_included)))#

filename_ex <- glue("/mnt/arcus/lab/users/kafadare/GAB_cohort/excluded/procedures_exclude_procDx_part{n}.rds")
print(paste0("Writing file: ", filename_ex))
saveRDS(procedures_excluded, filename_ex)

filename_inc <- glue("/mnt/arcus/lab/users/kafadare/GAB_cohort/included/procedures_include_procDx_part{n}.rds")
print(paste0("Writing file: ", filename_inc))
saveRDS(procedures_included, filename_inc)

if(upload == TRUE){
#Upload phecodes_excluded to SQL arcus
#procedures_excluded <- read.csv(filename)
bq_auth()
#proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]

proc_ord_ids_excluded <- procedures_excluded %>% select(proc_ord_id)

print(paste0("writing into SQL table now."))

bq_con <- dbConnect(bigrquery::bigquery(), project = proj_id)
dest_table = bq_table(project = proj_id, dataset = "lab", table = "proc_ord_id_exclude_procDx_EK")
bq_table_upload(
  x=dest_table, 
  values= proc_ord_ids_excluded, 
  create_disposition='CREATE_IF_NEEDED', 
  write_disposition='WRITE_APPEND')
}

#Look at distribution of GAB/age before & after exclusion of exclusion phecodes by date dx < date MRI, buffer 30 days.
if(plot == TRUE){
  GAB_Age_plot <- 
    procDx %>% distinct(proc_ord_id, .keep_all = TRUE) %>%
    ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
    scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))
  
  #look at distribution plot again
  GAB_Age_excluded_plot <- 
    procDx %>% 
    filter(!(proc_ord_id %in% procedures_excluded)) %>%
    distinct(proc_ord_id, .keep_all = TRUE) %>%
    ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
    scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))
}
