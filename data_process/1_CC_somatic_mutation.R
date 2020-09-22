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
    fread(file, quote = "")[, .(CHROM, POS, REF, ALT, T.AD, N.AD)] %>% unique
})
files <- dir("../data/CC_annovar_output/02_annovar/", "gene.hg19_multianno", full.names = T)
annGeneTabs <- llply(files , function(file) {
    fread(files[2], sep = ",") %>% 
        rename(c("Chr" = "CHROM", "Start" = "POS", "Ref" = "REF", "Alt" = "ALT")) %>% unique %>% 
        setkey(CHROM, POS, REF, ALT)
})
files <- dir("../data/CC_annovar_output/02_annovar/", "local.hg19_multianno", full.names = T)
annLocalTabs <- llply(files , function(file) {
    res = fread(file, sep = ",") %>% 
        rename(c("Chr" = "CHROM", "Start" = "POS", "Ref" = "REF", "Alt" = "ALT")) %>% unique
    if (is.logical(res$ALT[1])) {
        res$ALT = rep("T", length(res$ALT))
    }
    setkeyv(res, c("CHROM", "POS", "REF", "ALT"))
    res
})


parseAD <- function(data){
    ## input CHROM POS REF ALT C.AD M.AD N.AD
    ## return CHROM, POS, REF, ALT, C.F, C.A, M.F, M.A, N.F, N.A
    T.depth <- data$T.AD %>% sub("\\.", "0,0", .) %>%strsplit(",")%>% ldply(as.numeric) %>% set_names(c("T.R", "T.A"))
    N.depth <- data$N.AD %>% sub("\\.", "0,0", .) %>%strsplit(",")%>% ldply(as.numeric) %>% set_names(c("N.R", "N.A"))
    data.table(data[, .(CHROM, POS, REF, ALT)], T.depth, N.depth, key = c("CHROM", "POS", "REF", "ALT"))
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
