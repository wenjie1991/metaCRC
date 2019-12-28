# input         ../data/01_prepare_data.allele_frequency.output.tsv
# output        ../data/pyclone/tsv/

library(data.table)
library(magrittr)
library(readr)
library(plyr)

raw = fread("../data/01_prepare_data.allele_frequency.output.tsv")

personList = raw$personID %>% unique
names(raw)

prepare_tsv = function(person_tab) {
	person_tab_primary = person_tab[, .(
		mutation_id = mutation_id,
		ref_counts = C.R,
		var_counts = C.A,
		normal_cn = 2,
		minor_cn = 0,
		major_cn = 2
		)]
	person_tab_metastasis = person_tab[, .(
		mutation_id = mutation_id,
		ref_counts = M.R,
		var_counts = M.A,
		normal_cn = 2,
		minor_cn = 0,
		major_cn = 2
		)]
	list(primary = person_tab_primary, metastasis = person_tab_metastasis)
}

for (p in personList) {
	person_tab = raw[personID == p]
	person_prepared_tab = prepare_tsv(person_tab)
	write_tsv(person_prepared_tab$primary, paste0("../data/pyclone/tsv/", p, "_primary.tsv"))
	write_tsv(person_prepared_tab$metastasis, paste0("../data/pyclone/tsv/", p, "_metastasis.tsv"))
}


