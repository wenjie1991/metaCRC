# input             ../data/phenotype.csv, 
#                   ../data/00_mean_coverage.tsv
#                   ..data/patient_metrics.csv       
# output            Figure5A, Figure5B, Figure5C

library(stringr)
library(ggplot2)
library(magrittr)
library(data.table)

library(RColorBrewer)
colors = brewer.pal(5, "Blues")[c(3, 5)]
colors[3:5] =  brewer.pal(5, "Greens")[c(3:5)]

cd = fread("../data/phenotype.csv")
cd[grepl("^20", pathoID), patientID := paste0("ID", substr(pathoID, 3, 9))]


coverage = fread("../data/00_mean_coverage.tsv")
unqualify_sample = coverage[mean < 50 & !grepl("sn$", sample_name), sample_name] %>% str_match(., "(^\\w+)_") %>% extract(, 2) %>% unique


dat = fread("../data/patient_metrics.csv")
dat_sub1 = dat[grepl("^(AMC)|(ID11)", patientID)][, subset := paste0("1.All(", .N, ")")][!is.na(KSD)]
dat_sub2 = dat_sub1[!(grepl("^ID", patientID) & patientID %in% cd[Both_treated != "Chemonaive", patientID])][, subset := paste0("2.Remove treated(", .N, ")")]
dat_sub3 = dat_sub2[!patientID %in% unqualify_sample][, subset := paste0("3.Remove low depth(", .N, ")")]

dat_for_plot = rbind(dat_sub1, dat_sub2, dat_sub3)

dat_for_plot$patientID %>% unique
ggplot(dat_for_plot) + aes(x = KSD) + geom_density(aes(fill = subset), alpha = 0.5)
