# input     ../data/matched_triols_annovar_output/03_genotype_field/
#           ../data/matched_triols_annovar_output/02_annovar/
#           ../data/phenotype.csv
# output    ../data/01_trios_somaticMutation_raw.txt
#           ../data/01_trios_somaticMutation.txt
#           ../data/01_tumor_purity.tsv

library(data.table)
library(magrittr)
library(stringr)
library(readr)
library(plyr)
library(parallel)

#+ read data, include=FALSE
files <- dir("../data/matched_triols_annovar_output/03_genotype_field/", full.names = T)
sampleLable = str_match(files, "(ID|EV-|AMC|P-)\\d+")[, 1]
i.include = c(1:106, 122:137) # Included target gene seq, remove the MSI
sampleLable[setdiff(1:length(files), i.include)]
# The exclude samples
#  "ID0825371" "ID1119999" "ID185"     "ID250"     "ID262"     "ID278"
#  "ID353"     "ID381"     "ID413"     "ID503"     "ID509"     "ID523"
#  "ID526"     "ID627"     "ID707"     "ID718"     "ID721"
cd = fread("../data/phenotype.csv")
sampleID_nolivermeta = cd[metastasis_site != "Liver" | MS == "MSI", pathoID]
i.include = i.include[! (sampleLable[i.include] %in% sampleID_nolivermeta)]

#' Include sample
sampleLable[i.include]
sampleLable[setdiff(1:137, i.include)]

#+ read in mutation, eval=F
files <- dir("../data/matched_triols_annovar_output/03_genotype_field/", full.names = T)[i.include]
person <- str_match(files, "(ID|AMC|EV-|P-)\\d+")[, 1]

mutTabs <- llply(files, function(f) {
    fread(f)[, .(CHROM, POS, REF, ALT, C.AD, M.AD, N.AD)] %>% unique
})

files <- dir("../data/matched_triols_annovar_output/02_annovar/", "gene.hg19_multianno", full.names = T)[i.include]
annGeneTabs <- llply(files , function(file) {
    fread(file, sep = ",") %>% 
        rename(c("Chr" = "CHROM", "Start" = "POS", "Ref" = "REF", "Alt" = "ALT")) %>% unique %>% 
        setkey(CHROM, POS, REF, ALT)
})
files <- dir("../data/matched_triols_annovar_output/02_annovar/", "local.hg19_multianno", full.names = T)[i.include]
annLocalTabs <- llply(files , function(file) {
    fread(file, sep = ",") %>% 
        rename(c("Chr" = "CHROM", "Start" = "POS", "Ref" = "REF", "Alt" = "ALT")) %>% unique %>% 
        setkey(CHROM, POS, REF, ALT)
})


#+ functions and theme, include=F
parseAD <- function(data){
    ## input CHROM POS REF ALT C.AD M.AD N.AD
    ## return CHROM, POS, REF, ALT, C.F, C.A, M.F, M.A, N.F, N.A
    C <- data$C.AD %>% sub("\\.", "0,0", .) %>%strsplit(",")%>% ldply(as.numeric) %>% set_names(c("C.R", "C.A"))
    M <- data$M.AD %>% sub("\\.", "0,0", .) %>%strsplit(",")%>% ldply(as.numeric) %>% set_names(c("M.R", "M.A"))
    N <- data$N.AD %>% sub("\\.", "0,0", .) %>%strsplit(",")%>% ldply(as.numeric) %>% set_names(c("N.R", "N.A"))
    data.table(data[, .(CHROM, POS, REF, ALT)], C, M, N, key = c("CHROM", "POS", "REF", "ALT"))
}

#' # Somatic Identify Mutation
#+ somatic mutation Call, eval=F
cutpointT <- 0.05
cutpointN <- 0.01
identifySomaticMut <- function(data, type) {
    ## input CHROM POS REF ALT C.F C.A M.F M.F N.A 
    ## return CHROM, POS, REF, ALT, C.F, C.A, M.F, M.A, N.F, N.A
    switch(type,
        tumor = 
            data[
                (N.A / (N.A + N.R + 1) <= cutpointN) 
                & ((M.A / (M.A + M.R + 1) > cutpointT) | (C.A / (C.A + C.R + 1) > cutpointT))
                & (M.A > 3 | C.A > 3)
                & (N.A + N.R > 5)
                ]
        ,
        primary = 
            data[
                (N.A / (N.A + N.R + 1) <= cutpointN) 
                & (C.A / (C.A + C.R + 1) > cutpointT)
                & (C.A > 3)
                & (N.A + N.R > 5)
                ]
        ,
        metastasis = 
            data[
                (N.A / (N.A + N.R + 1) <= cutpointN) 
                & (M.A / (M.A + M.R + 1) > cutpointT) 
                & (M.A > 3)
                & (N.A + N.R > 5)
                ]
        )
}

parsedTabs <- mclapply(mutTabs, parseAD, mc.cores = 5)
normalMutTabs <- rbindlist(parsedTabs)[N.A / (N.A + N.R) > 0.05]
tumorMutTabs <- mclapply(parsedTabs, identifySomaticMut, type = "tumor")
primaryMutTabs <- mclapply(parsedTabs, identifySomaticMut, type = "primary")
metastasisMutTabs <- mclapply(parsedTabs, identifySomaticMut, type = "metastasis")

#' # Metastasis specific mutation identify
# add the annotation to the mutation data
annTabs <- llply(1:length(tumorMutTabs), function(i) {
    if (nrow(tumorMutTabs[[i]]) != 0)  {
        annLocalTabs[[i]]$ALT %<>% as.character
        annGeneTabs[[i]]$ALT %<>% as.character
        tumorMutTabs[[i]]$ALT %<>% as.character
        data.table(merge(merge(annLocalTabs[[i]], annGeneTabs[[i]], by = c("CHROM", "POS", "REF", "ALT")), tumorMutTabs[[i]], by = c("CHROM", "POS", "REF", "ALT")), sample=person[i])
    }
}) %>% rbindlist


phenotype <- fread("../data/phenotype.csv")
phenotype[, personID := pathoID]
phenotype[!grepl("EV|P-", pathoID), personID := paste0("ID", sub("^20", "", pathoID))]
setkey(phenotype, personID)
setkey(annTabs, sample)
phenotype[, personID]
annTabs <- phenotype[annTabs]
annTabs[is.na(MS), MS := "MSS"]
geneName <- annTabs$Gene.refGene %>%  strsplit(",")
annTabs <- annTabs[rep(1:nrow(annTabs), unlist(lapply(geneName, length)))]
annTabs$Gene.refGene <- unlist(geneName)

# output the raw mutation result
write_tsv(annTabs, "../data/01_trios_somaticMutation_raw.txt")


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

annTabs3 <- annTabs[
    (!(ExonicFunc.refGene %in% c("synonymous SNV")))
    & (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
    & (as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
    & (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
    & (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
    ]
# 
## Purity
prepare_af = function(af_alt, af_ref, af_n) {
    af_mut = af_alt + af_ref
    ratio_n = af_mut / af_n * sum(af_n) / sum(af_mut)
    cnv_i = abs(log2(ratio_n)) < log2(1.2)
    deepth = af_mut > 20 & af_mut < 100
    somatic_mut = (af_alt / (af_ref + af_alt + 1) > 0.15) & (af_alt / (af_ref + af_alt + 1) < 0.7);

    af_ref_s = af_ref[cnv_i & deepth & somatic_mut]
    af_alt_s = af_alt[cnv_i & deepth & somatic_mut]
    list(af_ref_s = af_ref_s, af_alt_s = af_alt_s)
}

personID_list = annTabs3$personID %>% unique
local_personID = personID_list %>% grep("^ID|AM", ., value=T)
 
person = local_personID[i]
library(parallel)
library(PurBayes)
library(plyr)
purity_l = mclapply(local_personID, function(person) {
    af_alt_p = annTabs3[personID == person, C.A]
    af_ref_p = annTabs3[personID == person, C.R]
    af_alt_m = annTabs3[personID == person, M.A]
    af_ref_m = annTabs3[personID == person, M.R]
    af_n     = annTabs3[personID == person, N.A + N.R]

    af_l_p = prepare_af(af_alt_p, af_ref_p, af_n)
    af_l_m = prepare_af(af_alt_m, af_ref_m, af_n)

    PB_p = try(PurBayes(af_l_p$af_ref_s + af_l_p$af_alt_s, af_l_p$af_alt_s), TRUE)
    PB_summary_p = summary(PB_p)
    PB_m = try(PurBayes(af_l_m$af_ref_s + af_l_m$af_alt_s, af_l_m$af_alt_s), TRUE)
    PB_summary_m = summary(PB_m)
    r = c(person, PB_summary_p[[1]], PB_summary_m[[1]])
    r
        }, mc.cores = 7)
purity_m = ldply(purity_l[-c(1, 13, 14)])
names(purity_m) = c("personID", "median_p", "p5_p", "p95_p", "median_m", "p5_m", "p95_m") 
write_tsv(purity_m, "../data/01_tumor_purity.tsv")

# Output filtered mutation
setkeyv(annTabs2, c('CHROM', 'POS', 'REF', 'ALT'))
setkeyv(normalMutTabs, c('CHROM', 'POS', 'REF', 'ALT'))
annTabs2 <- annTabs2[!normalMutTabs]
write_tsv(annTabs2, "../data/01_trios_somaticMutation.txt")
write_tsv(annTabs2[C.A / (C.A + C.R + 1) > 0.1 | M.A / (M.A + M.R + 1) > 0.1], "../data/01_trios_somaticMutation_upload.txt")
