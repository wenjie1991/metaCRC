# input         ../data/01_trios_somaticMutation_raw.txt
#               ../data/01_trios_somaticMutation.txt
#               ../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv
# output
#               ../data/02_VAF_t_test.tsv
#               ../data/02_Primary_mut_rate.tsv
#               ../data/02_Metastasis_mut_rate.tsv
#               ../data/02_Primary_mut_rate_right.tsv
#               ../data/02_Primary_mut_rate_left.tsv
#               ../data/02_Metastasis_mut_rate_right.tsv
#               ../data/02_Metastasis_mut_rate_left.tsv
#               ../data/02_Primary_mut_rate_syn.tsv
#               ../data/02_Primary_mut_rate_sub.tsv
#               ../data/02_Metastasis_mut_rate_syn.tsv
#               ../data/02_Metastasis_mut_rate_sub.tsv
#               ../data/02_Primary_mut_rate_both.tsv
#               ../data/02_Primary_mut_rate_naive.tsv
#               ../data/02_Primary_mut_rate_meta.tsv
#               ../data/02_Metastasis_mut_rate_both.tsv
#               ../data/02_Metastasis_mut_rate_naive.tsv
#               ../data/02_Metastasis_mut_rate_meta.tsv
#               ../data/02_Primary_mut_rate_MSK.tsv
#               ../data/02_Primary_mut_rate_Local.tsv
#               ../data/02_Primary_mut_rate_Oncotarget.tsv
#               ../data/02_Metastasis_mut_rate_MSK.tsv
#               ../data/02_Metastasis_mut_rate_Local.tsv
#               ../data/02_Metastasis_mut_rate_Oncotarget.tsv
#               ../data/02_mutDT.tsv

library(magrittr)
library(data.table)
library(plyr)
library(readr)

tumorVAFThreshold = 0.10
# normalVAFThreshold = 0.01 // the filter was applied in previous step

#+ Read in data
annTabs_raw <- fread("../data/01_trios_somaticMutation_raw.txt")
annTabs2.0 <- fread("../data/01_trios_somaticMutation.txt")
msk_impact_341 = readLines("../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv")
annTabs2.1 <- annTabs2.0[(metastasis_site %in% c("Liver")) | grepl("AMC", personID)]
annTabs2.2 <- annTabs2.1[MS == "MSS"]
annTabs2 = annTabs2.2[Gene.refGene %in% msk_impact_341]
n_metastasis = annTabs_raw$personID %>% table %>% length
n_right_metastasis = annTabs_raw[primary_site == "right"]$personID %>% table %>% length
n_left_metastasis = annTabs_raw[primary_site == "left"]$personID %>% table %>% length
n_syn_metastasis = annTabs_raw[Resection_timing == "Concurrent"]$personID %>% table %>% length
n_sub_metastasis = annTabs_raw[Resection_timing == "Subsequent"]$personID %>% table %>% length
n_naive_metastasis = annTabs_raw[Both_treated == "Chemonaive"]$personID %>% table %>% length
n_both_metastasis = annTabs_raw[Both_treated == "Both_treated"]$personID %>% table %>% length
n_meta_metastasis = annTabs_raw[Both_treated == "Metastasis_treated"]$personID %>% table %>% length
print(n_metastasis)

#' # Count the sample size
(annTabs_raw[, .(personID, primary_site)] %>% unique)[, primary_site] %>% table(useNA = "ifany")
(annTabs_raw[, .(personID, primary_site)] %>% unique)[, primary_site] %>% table(useNA = "ifany") %>% sum
(annTabs_raw[, .(personID, Both_treated)] %>% unique)[, Both_treated] %>% table(useNA = "ifany")
(annTabs_raw[, .(personID, Resection_timing)] %>% unique)[, Resection_timing] %>% table(useNA = "ifany")

annTabs_raw[, unique(personID) %>% sub("\\d+", "", .) %>% table]
annTabs_raw[, unique(personID)]



#' # Mutation count
#' ## Somatic mutation
somaticMutation <- annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) | (M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) 
    , .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]
(somaticMutation[n >= 4][, p := n / length(unique(annTabs2[, personID]))])

#' ## Primary Mutation
somaticMutationP <- annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2), .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]
(somaticMutationP[n >= 4])


#' ## Metastasis Mutation
somaticMutationM <- annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2), .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]
(somaticMutationM[n >= 4])

somaticMutationP[, group := "primary"]
somaticMutationM[, group := "metastasis"]
somaticMutationPM <- rbind(somaticMutationP, somaticMutationM)
somaticMutationPM$Gene.refGene %<>% factor(
    levels = somaticMutationPM[, .(s=sum(n)), by=Gene.refGene][order(s, decreasing=T), Gene.refGene])
somaticMutation$Gene.refGene %<>% factor(levels=somaticMutation$Gene.refGene)

#' ## Visualization
#' ### Mutation Candidate between stage
geneName <- somaticMutation[n >= 4, Gene.refGene]
AF.candidate <- ldply(as.character(geneName), function(i) {
    print(i)
    x <- annTabs2[Gene.refGene == i, .(AF.C = C.A / (C.A + C.R + 1) , AF.M = M.A / (M.A + M.R + 1))][, AF.M - AF.C]
    r = t.test(x)
    c(deltaMedian = mean(x),  p = r$p.value, symbol=i)
}) %>% data.table
AF.candidate[, p.adj := p.adjust(p, method="BH")]
(AF.candidate[p < 0.10])
candidateGene = somaticMutationPM[n > 4, Gene.refGene] %>% unique
write_tsv(AF.candidate, "../data/02_VAF_t_test.tsv")

exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_metastasis * 100
exonFreqC = exonFreq
write_tsv(exonFreqC, "../data/02_Primary_mut_rate.tsv")

exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]


exonFreq$freq =  exonFreq$freq / n_metastasis * 100
exonFreqM = exonFreq
write_tsv(exonFreqM, "../data/02_Metastasis_mut_rate.tsv")

#' ### Mutation between location
exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) & primary_site == "right", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_right_metastasis * 100
exonFreqC_right = exonFreq
write_tsv(exonFreqC_right, "../data/02_Primary_mut_rate_right.tsv")

exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) & primary_site == "left", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_left_metastasis * 100
exonFreqC_left = exonFreq
write_tsv(exonFreqC_left, "../data/02_Primary_mut_rate_left.tsv")


exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) & primary_site == "right", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_right_metastasis * 100
exonFreqM_right = exonFreq
write_tsv(exonFreqM_right, "../data/02_Metastasis_mut_rate_right.tsv")

exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) & primary_site == "left", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_left_metastasis * 100
exonFreqM_left = exonFreq
write_tsv(exonFreqM_left, "../data/02_Metastasis_mut_rate_left.tsv")

#' ### Mutation between resection time
exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) & Resection_timing == "Concurrent", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_syn_metastasis * 100
exonFreqC_syn = exonFreq
write_tsv(exonFreqC_syn, "../data/02_Primary_mut_rate_syn.tsv")

exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) & Resection_timing == "Subsequent", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_sub_metastasis * 100
exonFreqC_sub = exonFreq
write_tsv(exonFreqC_sub, "../data/02_Primary_mut_rate_sub.tsv")


exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) & Resection_timing == "Concurrent", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_syn_metastasis * 100
exonFreqM_syn = exonFreq
write_tsv(exonFreqM_syn, "../data/02_Metastasis_mut_rate_syn.tsv")

exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) & Resection_timing == "Subsequent", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_sub_metastasis * 100
exonFreqM_sub = exonFreq
write_tsv(exonFreqM_sub, "../data/02_Metastasis_mut_rate_sub.tsv")

#' ### Mutation between Chemotherapy 
exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) & Both_treated == "Both_treated", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_both_metastasis * 100
exonFreqC_both = exonFreq
write_tsv(exonFreqC_both, "../data/02_Primary_mut_rate_both.tsv")

exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) & Both_treated == "Chemonaive", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_naive_metastasis * 100
exonFreqC_naive = exonFreq
write_tsv(exonFreqC_naive, "../data/02_Primary_mut_rate_naive.tsv")

exonFreq = annTabs2[(C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2) & Both_treated == "Metastasis_treated", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq =  exonFreq$freq / n_meta_metastasis * 100
exonFreqC_meta = exonFreq
write_tsv(exonFreqC_meta, "../data/02_Primary_mut_rate_meta.tsv")

exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) & Both_treated == "Both_treated", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_both_metastasis * 100
exonFreqM_both = exonFreq
write_tsv(exonFreqM_both, "../data/02_Metastasis_mut_rate_both.tsv")

exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) & Both_treated == "Chemonaive", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_naive_metastasis * 100
exonFreqM_naive = exonFreq
write_tsv(exonFreqM_naive, "../data/02_Metastasis_mut_rate_naive.tsv")

exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2) & Both_treated == "Metastasis_treated", .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_meta_metastasis * 100
exonFreqM_meta = exonFreq
write_tsv(exonFreqM_meta, "../data/02_Metastasis_mut_rate_meta.tsv")


# n_genomeBiology = annTabs_raw[grepl("^EV", personID)]$personID %>% table %>% length
n_MSK = annTabs_raw[grepl("^EV", personID) | grepl("^P\\-", personID)]$personID %>% table %>% length
n_local = annTabs_raw[grepl("^ID\\d", personID)]$personID %>% table %>% length
n_oncotarget = annTabs_raw[grepl("^AMC", personID)]$personID %>% table %>% length
# n_NG = annTabs_raw[grepl("^P\\-", personID)]$personID %>% table %>% length
(data.frame(n_MSK, n_local, n_oncotarget))

exonFreq = annTabs2[(grepl("^EV", personID) | grepl("^P\\-", personID)) & (C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_MSK * 100
exonFreqM_meta = exonFreq
write_tsv(exonFreqM_meta, "../data/02_Primary_mut_rate_MSK.tsv")

exonFreq = annTabs2[grepl("^ID\\d", personID) & (C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_local * 100
exonFreqM_meta = exonFreq
write_tsv(exonFreqM_meta, "../data/02_Primary_mut_rate_Local.tsv")

exonFreq = annTabs2[grepl("^AMC", personID) & (C.A / (C.A + C.R + 1) > tumorVAFThreshold & C.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_oncotarget * 100
exonFreqM_meta = exonFreq
write_tsv(exonFreqM_meta, "../data/02_Primary_mut_rate_Oncotarget.tsv")

exonFreq = annTabs2[(grepl("^EV", personID) | grepl("^P\\-", personID)) & (M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_MSK * 100
exonFreqM_meta = exonFreq
write_tsv(exonFreqM_meta, "../data/02_Metastasis_mut_rate_MSK.tsv")

exonFreq = annTabs2[grepl("^ID\\d", personID) & (M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_local * 100
exonFreqM_meta = exonFreq
write_tsv(exonFreqM_meta, "../data/02_Metastasis_mut_rate_Local.tsv")

exonFreq = annTabs2[grepl("^AMC", personID) & (M.A / (M.A + M.R + 1) > tumorVAFThreshold & M.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq = exonFreq[, .(symbol, freq_n)]
exonFreq$freq = exonFreq$freq / n_oncotarget * 100
exonFreqM_meta = exonFreq
write_tsv(exonFreqM_meta, "../data/02_Metastasis_mut_rate_Oncotarget.tsv")

annTabs2$datasetName %>% table

#' ### Allelle Frequency 
deltaAF_all <- annTabs2[, .(
    personID
    , loci = paste0(CHROM, ":", POS, REF, "/", ALT)
    , cosmic70 = cosmic70
    , AAChange.refGene = AAChange.refGene
    , source= personID %>% sub("AMC\\d+", "AMC1", .) %>%  gsub("\\d", "+", .)
    , AF.N = N.A / (N.A + N.R + 1) %>% as.numeric
    , AF.M = M.A / (M.A + M.R + 1) %>% as.numeric
    , AF.C = C.A / (C.A + C.R + 1) %>% as.numeric
    , Gene.refGene)] %>% data.table
AF_all = deltaAF_all
AF_all[, purity_p := max(AF.C[Gene.refGene == "APC"]), by=personID]
AF_all[, purity_m := max(AF.M[Gene.refGene == "APC"]), by=personID]
AF_all = AF_all[purity_p > 0 & purity_m > 0]
AF_all[, deltaAF := (AF.M) - (AF.C)]
deltaAF = AF_all[Gene.refGene %in% candidateGene][purity_p > 0.20 & purity_m > 0.20]
deltaAF$source %<>% revalue(c("ID+++++++"="Local", "AMC+"="Lim_2015", "EV-+++"="Brannon_2014", "P-+++++++"="Yaeger_2018"))

deltaAFlong <- melt(deltaAF, id.vars=c("personID", "deltaAF", "Gene.refGene", "source", "loci")
    , measure.vars=c("AF.M", "AF.C", "AF.N")
    , value.name="AF", variable="AF_Type")
deltaAFlong$AF_Type %<>% factor(levels = c("AF.N", "AF.C", "AF.M"))
deltaAF$Gene.refGene %<>% factor(levels = exonFreqM[order(freq, decreasing = T), symbol])
deltaAFlong = deltaAFlong[, .(AF = max(AF), deltaAF = deltaAF[which.max(AF)], loci=loci[which.max(AF)]), by=c("Gene.refGene", "source", "personID", "AF_Type")]

deltaAFlong_draw = deltaAFlong[, .(Gene.refGene, deltaAF, source)] %>% unique
deltaAFlong_draw$Gene.refGene %<>% factor(levels = deltaAFlong_draw[, .N, by=Gene.refGene][order(N, decreasing=T), Gene.refGene])
deltaAFlong_draw[, wilcox.test(deltaAF)$p.value, by=Gene.refGene]

#' ### Mutation Landscape
mutDT = deltaAF_all[, .(sample = personID, gene = Gene.refGene, AF1 = AF.C, AF2 = AF.M, loci, cosmic70, AAChange.refGene)]
write_tsv(mutDT, "../data/02_mutDT.tsv")
