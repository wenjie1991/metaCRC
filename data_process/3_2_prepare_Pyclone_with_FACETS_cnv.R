# Input: ../data/01_prepare_data.allele_frequency.output2.tsv
# Output: ../data/pyclone/tsv_FACETS_cnv

annTabs = fread("../data/01_trios_somaticMutation_raw.txt")

## Pyclone input with CNV infered by FACETS
annTabs4 <- annTabs[
    #     (!(ExonicFunc.refGene %in% c("synonymous SNV")))
    (C.A / (C.A + C.R + 1) > 0.1 | M.A / (M.A + M.R + 1) > 0.1)
    & (C.A + C.R > 20 & M.A + M.R > 20)
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

    af_alt_p = annTabs4[personID == person, C.A]
    af_ref_p = annTabs4[personID == person, C.R]
    af_alt_m = annTabs4[personID == person, M.A]
    af_ref_m = annTabs4[personID == person, M.R]
    af_n     = annTabs4[personID == person, N.A + N.R]

    chr      = annTabs4[personID == person, CHROM] %>% sub("chr", "", .)
    loc      = annTabs4[personID == person, POS] 
    gene_name = annTabs4[personID == person, Gene.refGene]

    cnf_p    = fread(str_glue("../data/FACETS/cncf_{person}_primary.tsv"))
    cnf_m    = fread(str_glue("../data/FACETS/cncf_{person}_metastasis.tsv"))

    #TBD
    cnv_p = lapply(1:length(loc), function(i) {
        cnf_p[chrom == chr[i] & start < loc[i] & loc[i] < end, .(mutation_id = paste0(chr[i], "_", loc[i], "_", gene_name[i]), chr = chr[i], start = loc[i], tcn.em,  lcn.em, ref_counts = af_ref_p[i], var_counts = af_alt_p[i])]
        }) %>% ldply %>% data.table
    cnv_m = lapply(1:length(loc), function(i) {
        cnf_m[chrom == chr[i] & start < loc[i] & loc[i] < end, .(mutation_id = paste0(chr[i], "_", loc[i], "_", gene_name[i]), chr = chr[i], start = loc[i], tcn.em, lcn.em, ref_counts = af_ref_m[i], var_counts = af_alt_m[i])]
        }) %>% ldply %>% data.table

	person_tab_primary = cnv_p[, .(
		mutation_id,
		ref_counts,
		var_counts,
		normal_cn = 2,
		minor_cn = lcn.em,
		major_cn = tcn.em - lcn.em 
		)]

	person_tab_metastasis = cnv_m[, .(
		mutation_id,
		ref_counts,
		var_counts,
		normal_cn = 2,
		minor_cn = lcn.em,
		major_cn = tcn.em - lcn.em
		)]
	out = list(primary = person_tab_primary, metastasis = person_tab_metastasis)

	return(out)
}) -> annTabsX

names(annTabsX) = local_personID

sapply(annTabsX, function(x) {nrow(x)})

# output_file = "../data/01_prepare_data.allele_frequency_FACETCNV.output.tsv"

# ldply(annTabsX, .id = "personID") %>% write_tsv(output_file)

dir.create("../data/pyclone/tsv_FACETS_cnv/")
for (i in 1:length(local_personID)) {
    p = local_personID[i]
	write_tsv(annTabsX[[i]]$primary, paste0("../data/pyclone/tsv_FACETS_cnv/", p, "_primary.tsv"))
	write_tsv(annTabsX[[i]]$metastasis, paste0("../data/pyclone/tsv_FACETS_cnv/", p, "_metastasis.tsv"))
}

tumor_purity = fread("../data/tumor_purity_facet.tsv")

library(stringr)
for (i in 1:length(local_personID)) {
    p = local_personID[i]
    primary_purity = tumor_purity[personID == p, purity_p]
    metastasis_purity = tumor_purity[personID == p, purity_m]
    if (is.na(primary_purity) | is.na(metastasis_purity)) { next }
    cmd = str_glue("bash ./run_pyclone.sh {p} {primary_purity} {metastasis_purity}")
    print(cmd)
    system(cmd)
}

setkey(tumor_purity, "personID")

conf_dt = data.table(job_name = local_personID, tumor_purity[local_personID, .(purity_p, purity_m)])  %>% na.omit

write_tsv(conf_dt, "../data/pyclone/script/config.tsv")
