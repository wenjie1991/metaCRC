# input         ../data/01_trios_somaticMutation_raw.txt
#               ../data/01_trios_somaticMutation.txt")
# output        Figure2.A, Figure2.B
#               ../data/02_mutDT.tsv

library(magrittr)
library(data.table)
library(ggplot2)
library(plyr)
library(readr)

tumorPurityThreshold = 0.10
normalPurityThreshold = 0.01

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

#' # Count the sample size
(annTabs_raw[, .(personID, primary_site)] %>% unique)[, primary_site] %>% table(useNA = "ifany")
(annTabs_raw[, .(personID, primary_site)] %>% unique)[, primary_site] %>% table(useNA = "ifany") %>% sum
(annTabs_raw[, .(personID, Both_treated)] %>% unique)[, Both_treated] %>% table(useNA = "ifany")
(annTabs_raw[, .(personID, Resection_timing)] %>% unique)[, Resection_timing] %>% table(useNA = "ifany")
annTabs_raw[, unique(personID) %>% sub("\\d+", "", .) %>% table]
annTabs_raw[, unique(personID)]


somaticMutation <- annTabs2[(C.A / (C.A + C.R + 1) > tumorPurityThreshold & C.A > 2) | (M.A / (M.A + M.R + 1) > tumorPurityThreshold & M.A > 2) 
    , .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]

#' ## Primary Mutation
somaticMutationP <- annTabs2[(C.A / (C.A + C.R + 1) > tumorPurityThreshold & C.A > 2), .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]


somaticMutationM <- annTabs2[(M.A / (M.A + M.R + 1) > tumorPurityThreshold & M.A > 2), .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]

somaticMutationP[, group := "primary"]
somaticMutationM[, group := "metastasis"]
somaticMutationPM <- rbind(somaticMutationP, somaticMutationM)
somaticMutationPM$Gene.refGene %<>% factor(
    levels = somaticMutationPM[, .(s=sum(n)), by=Gene.refGene][order(s, decreasing=T), Gene.refGene])
somaticMutation$Gene.refGene %<>% factor(levels=somaticMutation$Gene.refGene)

#' # Candidate Gene
geneName <- somaticMutation[n >= 4, Gene.refGene]
AF.candidate <- ldply(as.character(geneName), function(i) {
    print(i)
    x <- annTabs2[Gene.refGene == i, .(AF.C = C.A / (C.A + C.R + 1) , AF.M = M.A / (M.A + M.R + 1))][, AF.M - AF.C]
    r = t.test(x)
    c(deltaMedian = mean(x),  p = r$p.value, symbol=i)
}) %>% data.table
AF.candidate[, p.adj := p.adjust(p, method="BH")]
candidateGene = somaticMutationPM[n > 4, Gene.refGene] %>% unique

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


exonFreq = annTabs2[(M.A / (M.A + M.R + 1) > tumorPurityThreshold & M.A > 2), .(symbol = Gene.refGene, freq_n = length(unique(personID))), by = Gene.refGene]
exonFreq$freq =  exonFreq$freq / n_metastasis * 100
exonFreqM = exonFreq

deltaAF$Gene.refGene %<>% factor(levels = exonFreqM[order(freq, decreasing = T), symbol])
deltaAFlong = deltaAFlong[, .(AF = max(AF), deltaAF = deltaAF[which.max(AF)], loci=loci[which.max(AF)]), by=c("Gene.refGene", "source", "personID", "AF_Type")]

#+ Allele Frequency Wilcox P lt 0.3,  fig.height=7, fig.width=13, dev='svg'
deltaAFlong_draw = deltaAFlong[, .(Gene.refGene, deltaAF, source)] %>% unique
deltaAFlong_draw$Gene.refGene %<>% factor(levels = deltaAFlong_draw[, .N, by=Gene.refGene][order(N, decreasing=T), Gene.refGene])
ggplot(data= deltaAFlong_draw %>% unique, aes(x=Gene.refGene, deltaAF)) + 
    geom_boxplot(alpha=0.2) + 
    geom_jitter(width=0.3, aes(color=source), alpha=0.3, size=3) + 
    xlab("Gene") + ylab("VAF_Metastasis - VAF_Primary")


#' ### Allele Frequency of small Wilcox P value
#+ AF by sample type,  fig.height=7, fig.width=10, dev='svg'
for (i in 1:length(candidateGene)){
    plot_AF_line(candidateGene[i]) %>% print
}
