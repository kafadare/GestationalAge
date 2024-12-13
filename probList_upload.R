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
rowcount_query <- 'SELECT count(*) FROM lab.proc_probList_EK'
rowcount <- bq_project_query(proj_id, rowcount_query) %>% bq_table_download()

n <- as.numeric(Sys.getenv (x = "SPLIT_ID"))
#n <- 1
print(paste0("SPLIT_ID is: ", n))


#read in the arguments passed on from the shell script
args <- commandArgs(trailingOnly = TRUE)

#Is argument length > 0?
if (length(args) > 0) {
  chunks <- as.numeric(args[1])
}else{
  n = 1
  chunks = 1
}

#chunks <- 50

print(paste0("NUMBER OF CHUNKS is: ", chunks))

size <- ceiling(rowcount/chunks)
print(paste0("SIZE is: ", size))

saveRDS <- FALSE

plot <- FALSE

offset <- n-1

filename <- glue("/mnt/arcus/lab/users/kafadare/GAB_cohort/excluded_procIDs/proc_ord_ids_exclude_probList_part{n}.csv")
print(paste0("Reading file: ", filename))

#Upload phecodes_excluded to SQL arcus
procedures_excluded <- read.csv(filename)
bq_auth()
proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]

proc_ord_ids_excluded <- procedures_excluded %>% select(proc_ord_id)

print(paste0("writing into SQL table now."))

bq_con <- dbConnect(bigrquery::bigquery(), project = proj_id)
dest_table = bq_table(project = proj_id, dataset = "lab", table = "proc_ord_id_exclude_probList_EK")
bq_table_upload(
  x=dest_table, 
  values= proc_ord_ids_excluded, 
  create_disposition='CREATE_IF_NEEDED', 
  write_disposition='WRITE_APPEND')

#Look at distribution of GAB/age before & after exclusion of exclusion phecodes by date dx < date MRI, buffer 30 days.
if(plot == TRUE){
  GAB_Age_plot <- 
    probList %>% distinct(proc_ord_id, .keep_all = TRUE) %>%
    ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
    scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))
  
  #look at distribution plot again
  GAB_Age_excluded_plot <- 
    probList %>% 
    filter(!(proc_ord_id %in% procedures_excluded)) %>%
    distinct(proc_ord_id, .keep_all = TRUE) %>%
    ggplot(., aes(x = age_at_scan_months, y = gestational_age_num)) + geom_bin2d(bins = 100)+
    scale_fill_viridis_c(option = "viridis", limits=(c(5, 2100)))
}
