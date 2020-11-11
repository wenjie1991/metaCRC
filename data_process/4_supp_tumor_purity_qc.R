library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

theme0 <- theme_bw() + theme(
    text = element_text(size = 15),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 60, vjust = 1, size = 15, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "right",
) 

# tumor purity estimated by Purebayes,
purebayes_purity = fread("../data/01_tumor_purity.tsv")

# APC VAF
mut_tb = fread("../data/01_trios_somaticMutation.txt")
apc_primary = mut_tb[Gene.refGene == "APC", .(personID, vaf_primary = C.A / (C.R + C.A))]
apc_metastasis = mut_tb[Gene.refGene == "APC", .(personID, vaf_metastasis = C.A / (M.R + M.A))]
apc_vaf = merge(apc_primary, apc_metastasis)


d = apc_vaf[vaf_primary > 0 & vaf_metastasis > 0 & vaf_primary < 0.5 & vaf_metastasis < 0.5, .(
    # d = apc_vaf[, .(
    personID  = c(personID %>% paste0("-p"), personID %>% paste0("-m")), 
    vaf = c(vaf_primary, vaf_metastasis)
    )] %>% unique
d_purity = purebayes_purity[, .(
    personID  = c(personID %>% paste0("-p"), personID %>% paste0("-m")), 
    median = c(median_p, median_m),
    p5 = c(p5_p, p5_m),
    p95 = c(p95_p, p95_m)
    )]
d = d[personID %in% d_purity$personID]
d = cbind(d, d_purity[match(d$personID, d_purity$personID), -1])

d_tmp = d[, .(vaf = max(vaf)), by = personID]
d$apc_vaf_max = F
setkeyv(d, c("personID", "vaf"))
d[d_tmp, apc_vaf_max := T]


g = ggplot(d) + aes(x = vaf * 2, y = median, color = apc_vaf_max)
p = g + geom_point() + 
    geom_errorbar(aes(ymin=p5, ymax=p95, width=0))
p + theme0 + xlim(0, 1) + ylim(0, 1)


# only include the normal copy number apc mutation
person_list = mut_tb$personID %>% unique
person_list = person_list[!grepl("ID0", person_list)]

lapply(person_list, function(person) {
    af_alt_p = mut_tb[personID == person, C.A]
    af_ref_p = mut_tb[personID == person, C.R]
    af_alt_m = mut_tb[personID == person, M.A]
    af_ref_m = mut_tb[personID == person, M.R]
    af_n     = mut_tb[personID == person, N.A + N.R]
    gene_list = mut_tb[personID == person, Gene.refGene]

    normal_cnv_gene = function(af_alt, af_ref, af_n, gene_list) {
        af_mut = af_alt + af_ref
        ratio_n = af_mut / af_n * sum(af_n) / sum(af_mut)
        cnv_i = abs(log2(ratio_n)) < log2(1.2)
        deepth = af_mut > 20 & af_mut < 100
        somatic_mut = (af_alt / (af_ref + af_alt + 1) > 0.15) & (af_alt / (af_ref + af_alt + 1) < 0.7);
        is_pass = cnv_i & deepth & somatic_mut

        data.table(gene = gene_list[is_pass], vaf = af_alt[is_pass] / (af_alt[is_pass] + af_ref[is_pass]))[gene == "APC", vaf]
    }

    mut_vaf_p = normal_cnv_gene(af_alt_p, af_ref_p, af_n, gene_list)
    mut_vaf_m = normal_cnv_gene(af_alt_m, af_ref_m, af_n, gene_list)

    data.table(
        personID = c(rep(paste0(person, "-p"), length(mut_vaf_p)), rep(paste0(person, "-m"), length(mut_vaf_m))),
        vaf = c(mut_vaf_p, mut_vaf_m)
        )
}) -> filtered_vaf

d_sub = filtered_vaf %>% rbindlist
d_sub = cbind(d_sub, d_purity[match(d_sub$personID, d_purity$personID), .(median, p5, p95)])

g = ggplot(d_sub) + aes(x = vaf, y = median)
p = g + geom_point() + geom_errorbar(aes(ymin=p5, ymax=p95, width=0))
p + theme0 


