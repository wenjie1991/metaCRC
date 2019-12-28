# Input         ../data/CC_annovar_output/03_genotype_field/
#               ../data/CC_annovar_output/02_annovar/
#               ../data/CC_annovar_output/02_annovar/
# Output        ../data/01_CC_somaticMutation_raw.txt
#               ../data/01_CC_somaticMutation.txt


library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(plyr)
library(parallel)
library(readr)

#' This script identify somatic mutations and annote them.
#' Include sample
files <- dir("../data/CC_annovar_output/03_genotype_field/", full.names = T)
person <- str_match(files, "(ID|AMC|EV-|P-)\\d+")[, 1]
mutTabs <- llply(files, function(file) {
    fread(file, sep = "\t")[, .(CHROM, POS, REF, ALT, T.AD, N.AD)] %>% unique
})
files <- dir("../data/CC_annovar_output/02_annovar/", "gene.hg19_multianno", full.names = T)
annGeneTabs <- llply(files , function(file) {
    fread(file, sep = ",") %>% 
        rename(c("Chr" = "CHROM", "Start" = "POS", "Ref" = "REF", "Alt" = "ALT")) %>% unique %>% 
        setkey(CHROM, POS, REF, ALT)
})
files <- dir("../data/CC_annovar_output/02_annovar/", "local.hg19_multianno", full.names = T)
annLocalTabs <- llply(files , function(file) {
    fread(file, sep = ",") %>% 
        rename(c("Chr" = "CHROM", "Start" = "POS", "Ref" = "REF", "Alt" = "ALT")) %>% unique %>% 
        setkey(CHROM, POS, REF, ALT)
})


parseAD <- function(data){
    ## input CHROM POS REF ALT C.AD M.AD N.AD
    ## return CHROM, POS, REF, ALT, C.F, C.A, M.F, M.A, N.F, N.A
    T <- data$T.AD %>% sub("\\.", "0,0", .) %>%strsplit(",")%>% ldply(as.numeric) %>% set_names(c("T.R", "T.A"))
    N <- data$N.AD %>% sub("\\.", "0,0", .) %>%strsplit(",")%>% ldply(as.numeric) %>% set_names(c("N.R", "N.A"))
    data.table(data[, .(CHROM, POS, REF, ALT)], T, N, key = c("CHROM", "POS", "REF", "ALT"))
}

#' # Somatic Identify Mutation
cutpointT <- 0.05
cutpointN <- 0.01
identifySomaticMut <- function(data, type) {
    ## input CHROM POS REF ALT C.F C.A M.F M.F N.A 
    ## return CHROM, POS, REF, ALT, C.F, C.A, M.F, M.A, N.F, N.A
    switch(type,
        tumor = 
            data[
                (N.A / (N.A + N.R + 1) <= cutpointN) 
                & (T.A > 3)
                & (N.A + N.R > 5)
                ]
        )
}

parsedTabs <- mclapply(mutTabs, parseAD, mc.cores = 5)
normalMutTabs <- rbindlist(parsedTabs)[N.A / (N.A + N.R + 1) > 0.05]
tumorMutTabs <- mclapply(parsedTabs, identifySomaticMut, type = "tumor")


#' # Metastasis specific mutation identify
# add the annotation to the mutation data
annTabs <- llply(1:length(tumorMutTabs), function(i) {
    data.table(annLocalTabs[[i]][annGeneTabs[[i]][tumorMutTabs[[i]]]], sample=person[i])
}) %>% rbindlist
phenotype <- fread("../data/CC_annovar_output/data_clinical_sample.txt", skip=4)
phenotype[, personID := PATIENT_ID]
phenotype[, MS := MSI_STATUS]
setkey(phenotype, personID)
setkey(annTabs, sample)
phenotype[, personID]
annTabs <- phenotype[annTabs]
geneName <- annTabs$Gene.refGene %>%  strsplit(",")
annTabs <- annTabs[rep(1:nrow(annTabs), unlist(lapply(geneName, length)))]
annTabs$Gene.refGene <- unlist(geneName)

# output the raw mutation result
write_tsv(annTabs, "../data/01_CC_somaticMutation_raw.txt")

# Filter the mutation by the effects of the mutations
annTabs2 <- annTabs[
    (!(ExonicFunc.refGene %in% c("synonymous SNV")))
    & (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
    & (as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
    & (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
    & (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
    & ((ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion", "stopgain"))
        | (SIFT_pred == "D" | Polyphen2_HDIV_pred == "D")
        | (Func.refGene %in% "splicing"))
    ]
 
# Output filtered mutation
setkeyv(annTabs2, c('CHROM', 'POS', 'REF', 'ALT'))
setkeyv(normalMutTabs, c('CHROM', 'POS', 'REF', 'ALT'))
annTabs2 <- annTabs2[!normalMutTabs]
write_tsv(annTabs2, "../data/01_CC_somaticMutation.txt")

