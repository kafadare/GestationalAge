require(bigrquery)
require(tidyr)
require(dplyr)

bq_auth()
proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]

my_query <- 'SELECT * FROM lab.GA_table_EK'
results <- bq_project_query(proj_id, my_query)
GA_dx_narratives <- bq_table_download(results, page_size=2000)

my_query <- 'SELECT * FROM lab.patient_phecode_dx'
results <- bq_project_query(proj_id, my_query)
patientPhecode_df <- bq_table_download(results, page_size=2000)

my_query <- 'SELECT * FROM arcus.patient'
results <- bq_project_query(proj_id, my_query)
patient_df <- bq_table_download(results, page_size=2000)

my_query <- 'SELECT * FROM arcus.procedure_order_narrative'
results <- bq_project_query(proj_id, my_query)
narratives_df <- bq_table_download(results, page_size=2000)

my_query <- 'SELECT * FROM arcus.procedure_order_diagnosis'
results <- bq_project_query(proj_id, my_query)
proc_ord_dx <- bq_table_download(results, page_size=2000)

my_query <- 'SELECT * FROM lab.grader_table_with_metadata'
results <- bq_project_query(proj_id, my_query)
graded.or.queued.includescoarsesearch <- bq_table_download(results, page_size=2000)

GAList <- patient_df[!is.na(patient_df$gestational_age),] #47591 patients

table(!(narratives_df$proc_ord_id %in% graded.or.queued.includescoarsesearch$proc_ord_id))
#shows how many narratives have not been evaluated yet
#note just focus on narratives_df for now but note that graders see impressions appended to narratives

#find reports not already graded
narratives_df_notgraded <- narratives_df[which(!narratives_df$proc_ord_id %in% 
                                                 graded.or.queued.includescoarsesearch$proc_ord_id),]

narratives_df_wDx <- merge(narratives_df_notgraded, proc_ord_dx, by = 'proc_ord_id')

GA_with_narratives <- merge(GAList, narratives_df_wDx, by = "pat_id")

##Diagnosis stuff now

my_query <- 'SELECT * FROM arcus.problem_list'
results <- bq_project_query(proj_id, my_query)
probList_df <- bq_table_download(results, page_size=2000)

#This takes forever to download
my_query <- 'SELECT * FROM arcus.encounter_diagnosis'
results <- bq_project_query(proj_id, my_query)
encounterDx_df <- bq_table_download(results, page_size=2000)

my_query <- 'SELECT * FROM arcus.master_diagnosis'
results <- bq_project_query(proj_id, my_query)
master_dx <- bq_table_download(results, page_size=2000)

###merge encounters & problem list with patient IDs from GAlist
full_table <- merge(GAList, encounterDx_df, by = "pat_id")
full_table <- merge(full_table, probList_df, by = "pat_id")


#encounterdx <- bq_table_download(results)
encounterdx <- encounterDx_df[,c('pat_id','dx_id')]
encounterdx2 <- paste0(encounterdx$pat_id,sep = '.isep.',encounterdx$dx_id)
encounter.dx.nodups <- encounterdx[which(!duplicated(encounterdx2)),]
rm(encounterdx2)
#encounter.dx.nodups is a matrix of pat_ids and dx_ids codes without duplicates

problistdx <- bq_table_download(results)
problistdx <- problistdx[,c('pat_id','dx_id')]
problistdx2 <- paste0(problistdx$pat_id,sep = '.isep.',problistdx$dx_id)
problist.dx.nodups <- problistdx[which(!duplicated(problistdx2)),]
rm(problistdx2)

#same thing for problem list

placeholderdx <- data.frame(pat_id=session_request_2024_02_with_metadata.filtered$pat_id, dx_id=2475657929)

combineddx <- rbind(encounter.dx.nodups,problist.dx.nodups,placeholderdx)
combineddx2 <- paste0(combineddx$pat_id,sep = '.isep.',combineddx$dx_id)
combined.dx.nodubs <- combineddx[which(!duplicated(combineddx2)),]
rm(combineddx2)
#same thing for combined

####code to phecode etc
#pathtoICDcode '/home/alexanderba/arcus-basics/R/ICDexclusion_code_March2024/'
#multiple ICD codes for each phecode
icd10tophecode <- read.csv('/home/alexanderba/arcus-basics/R/ICDexclusion_code_March2024/Phecode_map_v1_2_icd10cm_beta.csv') #downloaded Oct12 2023  from https://phewascatalog.org/phecodes_icd10
phecode_defs <- read.csv('/home/alexanderba/arcus-basics/R/ICDexclusion_code_March2024/phecode_definitions1.2.csv')
phecodes_exclusion <- read.csv('/home/alexanderba/arcus/shared/filter-by-dx/dx-filters/phecodes_with_exclusion_TS_and_AAB_19April2024.csv')
excluded_phecodes <- phecodes_exclusion[which(phecodes_exclusion$exclude_or_include_AAB_TS == 'exclude'),]
icdDx_toexclude <- icd10tophecode[which(icd10tophecode$phecode %in% excluded_phecodes$phecode),]
exclusionList <- master_dx$dx_id[which(master_dx$icd10_list %in% icdDx_toexclude$icd10cm)]





#narratives_df_notgraded_withdx <- merge(narratives_df_notgraded, proc_ord_dx, by = 'proc_ord_id')
#narratives_df_notgraded_withdx_macro <- narratives_df_notgraded_withdx[narratives_df_notgraded_withdx$dx_id %in% macrocephalylist,]









##Function to expand ranges to sequences
# range_seq <- function(range_str, step_size = 0.01) {
#   ranges <- strsplit(range_str, ",")[[1]]
#   sequences <- lapply(ranges, function(x) {
#     range_parts <- as.integer(strsplit(x, "-")[[1]])
#     seq(range_parts[1], range_parts[2], by = step_size)
#   })
#   unlist(sequences)
# }
# full_phecodes_excluded <- phecodes_exclusion$phecode_exclude_range
# #if nothing in range, take the single value from the phecode column
# full_phecodes_excluded[which(full_phecodes_excluded == "")] <- phecodes_exclusion[which(full_phecodes_excluded == ""),]$phecode
# excluded_withRange <- full_phecodes_excluded[grep("-", full_phecodes_excluded)]
# excluded_withOutRange <- full_phecodes_excluded[grep("-", full_phecodes_excluded, invert = TRUE)]
# test = sapply(excluded_withRange, function(x) range_seq(x[1])) %>% unlist(., recursive = TRUE, use.names = FALSE) %>% unique(.)
# full_phecodes_excluded <- c(test, excluded_withOutRange) %>% unique(.)

#Will change this to exclude based on dx, instead of include based on dx
## get encounter where dx was in macrocephaly list above
## note that found the dx_ids by hand
my_query <- 'SELECT * 
FROM arcus.encounter_diagnosis
INNER JOIN arcus.encounter
ON encounter.encounter_id = encounter_diagnosis.encounter_id
WHERE dx_id IN ("3717864987","331880719","453026982","1265559707","100011086","1685423373","660163900","1110870303","1625765175","3018752883","491946354","2557658056","1386182681","1638349569","2633083054","2076433868","2137984412","1603475265","4162158048","3818776794","2924634774","3072552928","2835437960","970825525","913824043","2997260075","4256034044","3519859477","1510184");'
results <- bq_project_query(proj_id, my_query)
encounter.macro <- bq_table_download(results)

##get a proc table with encounter info
my_query <- 'SELECT *
FROM arcus.procedure_order
INNER JOIN arcus.encounter
ON procedure_order.encounter_id = encounter.encounter_id;'
results <- bq_project_query(proj_id, my_query)
proc_ord_age <- bq_table_download(results)

narratives_df_with_patid <- merge(narratives_df_notgraded, proc_ord_age)

encounter.macro.withMRI <- encounter.macro[which(encounter.macro$pat_id %in% narratives_df_with_patid$pat_id),]
narratives_df_with_patid_with_macro <-  narratives_df_with_patid[which(narratives_df_with_patid$pat_id %in% encounter.macro$pat_id),]

index.with.macro <- vector()
for(i in 1:nrow(narratives_df_with_patid_with_macro)){
  encounter.macro.withMR.pat_i <- encounter.macro.withMRI[encounter.macro.withMRI$pat_id== narratives_df_with_patid_with_macro$pat_id[i],]
  if(any(abs(encounter.macro.withMR.pat_i$effective_age-narratives_df_with_patid_with_macro$effective_age[i]) < 91)){
    index.with.macro <- c(index.with.macro, i) 
  }
 if(!i %% 200){ cat(i, ' ')}
}

narratives_df_with_patid_with_macro_90days <- narratives_df_with_patid_with_macro[index.with.macro,] #note this list subsumes above list
#this is all patients have a dx of marocephaly within 3 months of their scan