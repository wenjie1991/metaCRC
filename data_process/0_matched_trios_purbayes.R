# Input         ../data/01_trios_somaticMutation_raw.txt
# Output        ../data/01_tumor_purity.tsv
#               ../data/01_tumor_purity_syn.tsv
#               ../data/01_tumor_purity_facetscnv.tsv
#               ../data/01_tumor_purity_facetscnvsyn.tsv

annTabs = fread("../data/01_trios_somaticMutation_raw.txt")

annTabs4 <- annTabs[
    #     (!(ExonicFunc.refGene %in% c("synonymous SNV")))
    C.A / (C.A + C.R + 1) > 0.1 | M.A / (M.A + M.R + 1) > 0.1
    & (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
    & (as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
    & (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
    & (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
    ]

annTabs3 <- annTabs[
    C.A / (C.A + C.R + 1) > 0.1 | M.A / (M.A + M.R + 1) > 0.1
    & (!(ExonicFunc.refGene %in% c("synonymous SNV")))
    & (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
    & (as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
    & (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
    & (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
    ]
 
## Purity legacy
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
i = sapply(purity_l, length) == 7
purity_m = ldply(purity_l[i])
names(purity_m) = c("personID", "median_p", "p5_p", "p95_p", "median_m", "p5_m", "p95_m") 
write_tsv(purity_m, "../data/01_tumor_purity.tsv")
 

## Purity infered by data include synonymous SNV
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

personID_list = annTabs4$personID %>% unique
local_personID = personID_list %>% grep("^ID|AM", ., value=T)
 
library(parallel)
library(PurBayes)
library(plyr)
purity_l = mclapply(local_personID, function(person) {
    af_alt_p = annTabs4[personID == person, C.A]
    af_ref_p = annTabs4[personID == person, C.R]
    af_alt_m = annTabs4[personID == person, M.A]
    af_ref_m = annTabs4[personID == person, M.R]
    af_n     = annTabs4[personID == person, N.A + N.R]

    af_l_p = prepare_af(af_alt_p, af_ref_p, af_n)
    af_l_m = prepare_af(af_alt_m, af_ref_m, af_n)

    PB_p = try(PurBayes(af_l_p$af_ref_s + af_l_p$af_alt_s, af_l_p$af_alt_s), TRUE)
    PB_summary_p = summary(PB_p)
    PB_m = try(PurBayes(af_l_m$af_ref_s + af_l_m$af_alt_s, af_l_m$af_alt_s), TRUE)
    PB_summary_m = summary(PB_m)
    r = c(person, PB_summary_p[[1]], PB_summary_m[[1]])
    r
        }, mc.cores = 7)
i = sapply(purity_l, length) == 7
purity_m = ldply(purity_l[i])
names(purity_m) = c("personID", "median_p", "p5_p", "p95_p", "median_m", "p5_m", "p95_m") 
write_tsv(purity_m, "../data/01_tumor_purity_syn.tsv")

## Purity with the CNV inform from FACETS
prepare_af = function(cnf, af_alt, af_ref, af_n, chr, loc) {
    cnv_i = sapply(1:length(af_alt), function(i) {
        cnf[chrom == chr[i] & start < loc[i] & loc[i] < end, tcn.em == 2 & lcn.em == 1]
        }) %>% unlist

    af_ref_s = af_ref[cnv_i]
    af_alt_s = af_alt[cnv_i]
    list(af_ref_s = af_ref_s, af_alt_s = af_alt_s)
}

personID_list = annTabs4$personID %>% unique
local_personID = personID_list %>% grep("^ID|AM", ., value=T)
 
library(parallel)
library(stringr)
library(PurBayes)
library(plyr)
purity_l = mclapply(local_personID, function(person) {
    af_alt_p = annTabs3[personID == person, C.A]
    af_ref_p = annTabs3[personID == person, C.R]
    af_alt_m = annTabs3[personID == person, M.A]
    af_ref_m = annTabs3[personID == person, M.R]
    af_n     = annTabs3[personID == person, N.A + N.R]

    chrom    = annTabs3[personID == person, CHROM] %>% sub("chr", "", .)
    loc      = annTabs3[personID == person, POS]

    cnf_p = fread(str_glue("../data/FACETS/cncf_{person}_primary.tsv"))
    cnf_m = fread(str_glue("../data/FACETS/cncf_{person}_metastasis.tsv"))

    af_l_p = prepare_af(cnf_p, af_alt_p, af_ref_p, af_n, chrom, loc)
    af_l_m = prepare_af(cnf_m, af_alt_m, af_ref_m, af_n, chrom, loc)

    PB_p = try(PurBayes(af_l_p$af_ref_s + af_l_p$af_alt_s, af_l_p$af_alt_s), TRUE)
    PB_summary_p = summary(PB_p)
    PB_m = try(PurBayes(af_l_m$af_ref_s + af_l_m$af_alt_s, af_l_m$af_alt_s), TRUE)
    PB_summary_m = summary(PB_m)
    r = c(person, PB_summary_p[[1]], PB_summary_m[[1]])
    r
        }, mc.cores = 7)
i = sapply(purity_l, length) == 7
purity_m = ldply(purity_l[i])
names(purity_m) = c("personID", "median_p", "p5_p", "p95_p", "median_m", "p5_m", "p95_m") 
write_tsv(purity_m, "../data/01_tumor_purity_facetscnv.tsv")


## Purity with the CNV inform from FACETS with synonymous mutation
prepare_af = function(cnf, af_alt, af_ref, af_n, chr, loc) {
    cnv_i = sapply(1:length(af_alt), function(i) {
        cnf[chrom == chr[i] & start < loc[i] & loc[i] < end, tcn.em == 2 & lcn.em == 1]
        }) %>% unlist

    af_ref_s = af_ref[cnv_i]
    af_alt_s = af_alt[cnv_i]
    list(af_ref_s = af_ref_s, af_alt_s = af_alt_s)
}

personID_list = annTabs4$personID %>% unique
local_personID = personID_list %>% grep("^ID|AM", ., value=T)
 
library(parallel)
library(stringr)
library(PurBayes)
library(plyr)
purity_l = mclapply(local_personID, function(person) {
    af_alt_p = annTabs4[personID == person, C.A]
    af_ref_p = annTabs4[personID == person, C.R]
    af_alt_m = annTabs4[personID == person, M.A]
    af_ref_m = annTabs4[personID == person, M.R]
    af_n     = annTabs4[personID == person, N.A + N.R]

    chrom    = annTabs4[personID == person, CHROM] %>% sub("chr", "", .)
    loc      = annTabs4[personID == person, POS]


    cnf_p = fread(str_glue("../data/FACETS/cncf_{person}_primary.tsv"))
    cnf_m = fread(str_glue("../data/FACETS/cncf_{person}_metastasis.tsv"))

    af_l_p = prepare_af(cnf_p, af_alt_p, af_ref_p, af_n, chrom, loc)
    af_l_m = prepare_af(cnf_m, af_alt_m, af_ref_m, af_n, chrom, loc)

    PB_p = try(PurBayes(af_l_p$af_ref_s + af_l_p$af_alt_s, af_l_p$af_alt_s), TRUE)
    PB_summary_p = summary(PB_p)
    PB_m = try(PurBayes(af_l_m$af_ref_s + af_l_m$af_alt_s, af_l_m$af_alt_s), TRUE)
    PB_summary_m = summary(PB_m)
    r = c(person, PB_summary_p[[1]], PB_summary_m[[1]])
    r
        }, mc.cores = 7)
i = sapply(purity_l, length) == 7
purity_m = ldply(purity_l[i])
names(purity_m) = c("personID", "median_p", "p5_p", "p95_p", "median_m", "p5_m", "p95_m") 
write_tsv(purity_m, "../data/01_tumor_purity_facetscnvsyn.tsv")

