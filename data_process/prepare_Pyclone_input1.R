# input         ../data/01_trios_somaticMutation_raw.txt
# output        ../data/01_prepare_data.allele_frequency.output.tsv

annTabs = fread("../data/01_trios_somaticMutation_raw.txt")

annTabs3 <- annTabs[
	#     (!(ExonicFunc.refGene %in% c("synonymous SNV")))
	#     & (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
	(as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
	& (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
	& (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
	]


tumor_purity = fread("../data/01_tumor_purity.tsv")

personID_list = annTabs3$personID %>% unique
local_personID = personID_list %>% grep("^ID|AM", ., value=T)

library(parallel)
library(plyr)

lapply(local_personID, function(person) {

    af_alt_p = annTabs3[personID == person, C.A]
    af_ref_p = annTabs3[personID == person, C.R]
    af_alt_m = annTabs3[personID == person, M.A]
    af_ref_m = annTabs3[personID == person, M.R]
    af_n     = annTabs3[personID == person, N.A + N.R]

	purity_p = tumor_purity[personID == person, median_p]
	purity_m = tumor_purity[personID == person, median_m]

	af_p = af_alt_p + af_ref_p 
	af_m = af_alt_m + af_ref_m 

	ratio_p = af_p / af_n * sum(af_n) / sum(af_p)
	ratio_m = af_m / af_n * sum(af_n) / sum(af_m)
	cnv_i = abs(log2(ratio_p)) < log2(1.2) & abs(log2(ratio_m)) < log2(1.2)
	#     deepth_i = af_p > 20 & af_p < 300 & af_m > 20 & af_m < 300
	i = cnv_i #& deepth_i

	out = annTabs3[personID == person][i][, .(mutation_id = paste0(Gene.refGene, "_", CHROM, ":", POS, "-", End, "_", REF, ":", ALT), C.A, C.R, M.A, M.R, N.A, N.R)]

	return(out)
}) -> annTabs4
names(annTabs4) = local_personID

sapply(annTabs4, function(x) {nrow(x)})

ldply(annTabs4, .id = "personID") %>% write_tsv("../data/01_prepare_data.allele_frequency.output.tsv")
