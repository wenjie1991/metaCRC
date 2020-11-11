# input         ../data/01_trios_somaticMutation_raw.txt
# output        ../data/01_prepare_data.allele_frequency.output2.tsv

annTabs = fread("../data/01_trios_somaticMutation_raw.txt")


## Pyclone input with CNV infered by FACETS
annTabs4 <- annTabs[
    #     (!(ExonicFunc.refGene %in% c("synonymous SNV")))
    C.A / (C.A + C.R + 1) > 0.1 | M.A / (M.A + M.R + 1) > 0.1
    & (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
    & (as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
    & (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
    & (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
    ]

personID_list = annTabs4$personID %>% unique
local_personID = personID_list %>% grep("^ID|AM", ., value=T)

library(parallel)
library(plyr)
library(stringr)


lapply(local_personID, function(person) {
    person = personID_list[1]

    af_alt_p = annTabs4[personID == person, C.A]
    af_ref_p = annTabs4[personID == person, C.R]
    af_alt_m = annTabs4[personID == person, M.A]
    af_ref_m = annTabs4[personID == person, M.R]
    af_n     = annTabs4[personID == person, N.A + N.R]

    chr    = annTabs4[personID == person, CHROM] %>% sub("chr", "", .)
    loc      = annTabs4$POS

    cnf_p = fread(str_glue("../data/FACETS/cncf_{person}_primary.tsv"))
    cnf_m = fread(str_glue("../data/FACETS/cncf_{person}_metastasis.tsv"))

    #TBD
    cnv_p_l = lapply(1:length(loc), function(i) {
        cnf_p[chrom == chr[i] & start < loc[i] & loc[i] < end, .(tcn.em,  lcn.em)]
        }) 
    cnv_m_i = lapply(1:length(loc), function(i) {
        cnf_m[chrom == chr[i] & start < loc[i] & loc[i] < end, .(tcn.em, lcn.em)]
        })

	out = annTabs4[personID == person][!is.na(cnv_p_l)][, .(
        mutation_id = paste0(Gene.refGene, "_", CHROM, ":", POS, "-", End.x, "_", REF, ":", ALT)
        , C.A
        , C.R
        , M.A
        , M.R
        , N.A
        , N.R
        , tcn_p = cnv_p$tcn.em
        , lcn_p = cnv_p$lcn.em
        , tcn_m = cnv_m$tcn.em
        , lcn_p = cnv_m$lcn.em
        )]

	return(out)
}) -> annTabsX

names(annTabsX) = local_personID

sapply(annTabsX, function(x) {nrow(x)})

ldply(annTabsX, .id = "personID") %>% write_tsv("../data/01_prepare_data.allele_frequency.output2.tsv")
