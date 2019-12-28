# input         ../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv
#               ../data/CC_annovar_output/data_clinical_sample.txt
#               ../data/01_CC_somaticMutation_raw.txt
#               ../data/01_CC_somaticMutation.txt
# output
#               ../data/02_early_stage_CC_mut_rate.tsv
#               ../data/02_meta_stage_CC_mut_rate.tsv
#               ../data/02_meta_CC_mut_rate.tsv

#+ setup, include=F
library(magrittr)
library(data.table)
library(readr)
library(plyr)

cutPoint_T = 0.10
cutPoint_N = 0.01

#+ Read in data
msk_impact_341 = readLines("../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv")

cd = fread("../data/CC_annovar_output/data_clinical_sample.txt", skip=4)
annTabs_all <- fread("../data/01_CC_somaticMutation_raw.txt")

annTabs = fread("../data/01_CC_somaticMutation.txt")
annTabs = annTabs[Gene.refGene %in% msk_impact_341]

is_single = function(x, y) {
    d = data.table(x, y) %>% unique
    tab = tapply(d$y, d$x, length)
    names(tab)[tab == 1]
}


# Early stage samples
annTabs_early_stage_all = annTabs_all[MSI_STATUS == "MSS" & SAMPLE_TYPE == "Primary" & STAGE_AT_DIAGNOSIS %in% c("I", "II")][PATIENT_ID %in% is_single(PATIENT_ID, SAMPLE_ID)]
n_early_stage = annTabs_early_stage_all[, length(unique(SAMPLE_ID))]
annTabs_early_stage = annTabs[MSI_STATUS == "MSS" & SAMPLE_TYPE == "Primary" & STAGE_AT_DIAGNOSIS %in% c("I", "II")][PATIENT_ID %in% is_single(PATIENT_ID, SAMPLE_ID)]
# Metastasis samples
annTabs_meta_stage_all = annTabs_all[MSI_STATUS == "MSS" & SAMPLE_TYPE == "Primary" & STAGE_AT_DIAGNOSIS %in% c("IV")][PATIENT_ID %in% is_single(PATIENT_ID, SAMPLE_ID)]
n_meta_stage = annTabs_meta_stage_all[, length(unique(SAMPLE_ID))]
annTabs_meta_stage = annTabs[MSI_STATUS == "MSS" & SAMPLE_TYPE == "Primary" & STAGE_AT_DIAGNOSIS %in% c("IV")][PATIENT_ID %in% is_single(PATIENT_ID, SAMPLE_ID)]
# Liver metastasis sample
annTabs_liver_all = annTabs_all[MSI_STATUS == "MSS" & SAMPLE_TYPE == "Metastasis" & METASTATIC_BIOPSY_SITE == "Liver"][PATIENT_ID %in% is_single(PATIENT_ID, SAMPLE_ID)]
n_liver = annTabs_liver_all[, length(unique(SAMPLE_ID))]
annTabs_liver = annTabs[MSI_STATUS == "MSS" & SAMPLE_TYPE == "Metastasis" & METASTATIC_BIOPSY_SITE == "Liver"][PATIENT_ID %in% is_single(PATIENT_ID, SAMPLE_ID)]

#' # Count sample size
(annTabs_early_stage_all[,.(personID, CHEMO_EXP_SEQ_SPECIMEN )] %>% unique)[, table(CHEMO_EXP_SEQ_SPECIMEN, useNA="ifany" )]
(annTabs_meta_stage_all[,.(personID, CHEMO_EXP_SEQ_SPECIMEN  )] %>% unique)[, table(CHEMO_EXP_SEQ_SPECIMEN, useNA="ifany" )]
(annTabs_liver_all[,.(personID, CHEMO_EXP_SEQ_SPECIMEN       )] %>% unique)[, table(CHEMO_EXP_SEQ_SPECIMEN, useNA="ifany" )]
(annTabs_early_stage_all[,.(personID, STAGE_AT_DIAGNOSIS     )] %>% unique)[, table(STAGE_AT_DIAGNOSIS, useNA="ifany"     )]
(annTabs_meta_stage_all[,.(personID, STAGE_AT_DIAGNOSIS      )] %>% unique)[, table(STAGE_AT_DIAGNOSIS, useNA="ifany"     )]
(annTabs_liver_all[,.(personID, STAGE_AT_DIAGNOSIS           )] %>% unique)[, table(STAGE_AT_DIAGNOSIS, useNA="ifany"     )]
(annTabs_early_stage_all[,.(personID, PRIMARY_TUMOR_LOCATION )] %>% unique)[, table(PRIMARY_TUMOR_LOCATION, useNA="ifany" )]
(annTabs_meta_stage_all[,.(personID, PRIMARY_TUMOR_LOCATION  )] %>% unique)[, table(PRIMARY_TUMOR_LOCATION, useNA="ifany" )]
(annTabs_liver_all[,.(personID, PRIMARY_TUMOR_LOCATION       )] %>% unique)[, table(PRIMARY_TUMOR_LOCATION, useNA="ifany" )]

#' # Mutation
get_mut_freq = function(n, annTabs, filename) {
    mutFreq = annTabs[T.A / (T.R + T.A) > cutPoint_T & N.A / (N.A + N.R) < cutPoint_N, 
        .(N = length(unique(personID))), by=Gene.refGene][, .(symbol=Gene.refGene, freq=N / n * 100, freq_n=N)]
    write_tsv(mutFreq, filename)
}

#' ## meta stage mutation
get_mut_freq(n_early_stage, annTabs_early_stage, "../data/02_early_stage_CC_mut_rate.tsv")

#' ## early stage mutation
get_mut_freq(n_meta_stage, annTabs_meta_stage, "../data/02_meta_stage_CC_mut_rate.tsv")

#' ## liver mutation
get_mut_freq(n_liver, annTabs_liver, "../data/02_meta_CC_mut_rate.tsv")
