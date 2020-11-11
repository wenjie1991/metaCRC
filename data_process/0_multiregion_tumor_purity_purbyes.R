# input         ../data/multiregion_biopsies/03_annotation/
# output        ../data/multiregion_purbyes_estimateion.tsv

library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(readr)
library(plyr)
library(parallel)
library(PurBayes)


# VAF files
files = dir("../data/multiregion_biopsies/03_annotation/", pattern="tsv", full=T)
person = str_match(basename(files), "\\w+")[, 1]

calc_tumor_purity = function(af_ref, af_alt) {

    af_n = af_alt + af_ref
    af_n_all = colSums(af_n)
    normal_i = grep("ND", names(af_n))
    ratio_n = af_n / af_n[, normal_i, with=F][[1]] * af_n_all[normal_i]
    ratio_n = data.table(t(t(ratio_n) / af_n_all))

    #     hist(log2(ratio_n[[3]][abs(log2(ratio_n[[3]])) < log2(1.2)]))

    normal_i = grep("ND", names(af_ref))
    germline_mut = af_alt[, normal_i, with=F] / (af_ref[, normal_i, with=F] + af_alt[, normal_i, with=F] + 1) > 0.05
    germline_mut %<>% as.vector
    somatic_mut = (af_alt / (af_ref + af_alt + 1) > 0.15) & (af_alt / (af_ref + af_alt + 1) < 0.7);
    deepth = (af_alt + af_ref) > 20 & (af_alt + af_ref) < 100
    
    purity = mclapply((1:ncol(af_ref))[-normal_i], function(i) {

        cnv_i = abs(log2(ratio_n[[i]])) < log2(1.2)

        af_ref_s = af_ref[!germline_mut & somatic_mut[, i] & deepth[, i] & cnv_i][[i]]
        af_alt_s = af_alt[!germline_mut & somatic_mut[, i] & deepth[, i] & cnv_i][[i]]
        #         hist(af_alt_s / (af_alt_s + af_ref_s + 1))
        PB = PurBayes(af_ref_s + af_alt_s, af_alt_s)
        PB_summary = summary(PB)
        PB_summary[[2]]
    }, mc.cores=7)
    purity
}


read_annoted_vaf = function(i) {

    d = fread(files[i])
    af = d[, -c(1:5), with=F]
    af_ref = af[, 1:(ncol(af)/2), with=F]
    af_alt = af[, -c(1:(ncol(af)/2)), with=F]
    vaf_frac = af_alt / (af_ref + af_alt + 1)
    names(vaf_frac) = str_match(names(vaf_frac), "\\w+")[, 1] %>% paste0(".FRAC")
    vaf = data.table(d[, c(1:5), with=F], vaf_frac, af)
    vaf = rename(vaf, c("CHROM" = "Chr", "POS_start" = "Start", "POS_end" = "End", "REF" = "Ref", "ALT" = "Alt"))
    setkeyv(vaf, c("Chr", "Start", "End", "Ref", "Alt"))


    local_ann = fread(paste0("../data/multiregion_biopsies/03_annotation/", person[i], "_local.hg19_multianno.csv"), key=c("Chr", "Start", "End", "Ref", "Alt"))
    gene_ann = fread(paste0("../data/multiregion_biopsies/03_annotation/", person[i], "_gene.hg19_multianno.csv"), key=c("Chr", "Start", "End", "Ref", "Alt"))
    ann = merge(local_ann, gene_ann)

    
    annTabs = ann
    annTabs2 <- annTabs[
        (!(ExonicFunc.refGene %in% c("synonymous SNV")))
        & (Func.refGene %in% c("exonic", "exonic;exonic", "exonic;splicing", "ncRNA_exonic", "splicing"))
        & (as.numeric(esp6500siv2_all) < 0.01 | esp6500siv2_all == "." | is.na(esp6500siv2_all)) 
        & (as.numeric(`1000g2015aug_all`) < 0.01 | `1000g2015aug_all` == "." | is.na(`1000g2015aug_all`))
        & (as.numeric(ExAC_ALL) < 0.05 | ExAC_ALL == "." | is.na(ExAC_ALL))
        #         & ((ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion", "stopgain"))
        #             | (SIFT_pred == "D" | Polyphen2_HDIV_pred == "D")
        #             | (Func.refGene %in% "splicing"))
        ]

    combind = vaf[annTabs2]
    vaf_sub = combind[, 1:(ncol(vaf)), with=F]
    ann_sub = combind[, -c(1:(ncol(vaf))), with=F][, .(Gene.refGene, Func.refGene, ExonicFunc.refGene, AAChange.refGene)]

    ## Purity
    #     chr_list = c("chr2", "chr3", "chr5", "chr6", "chr9", "chr10", "chr11", "chr12", "chr14", "chr16")
    af_ref_sub = vaf_sub[, grep("REF", names(vaf_sub)), with=F]
    af_alt_sub = vaf_sub[, grep("ALT", names(vaf_sub)), with=F]
    #     af_ref_sub = vaf[Chr %in% chr_list, grep("REF", names(vaf_sub)), with=F]
    #     af_alt_sub = vaf[Chr %in% chr_list, grep("ALT", names(vaf_sub)), with=F]
    #     af_ref_sub = vaf[, grep("REF", names(vaf_sub)), with=F]
    #     af_alt_sub = vaf[, grep("ALT", names(vaf_sub)), with=F]



    tumor_purity_l = calc_tumor_purity(af_ref=af_ref_sub, af_alt=af_alt_sub)

    vaf_frac_sub = vaf_sub[, c(1:5, grep("FRAC", names(vaf_sub))), with=F]
    m_sub = data.table(vaf_frac_sub, ann_sub)


    p1 = lapply(tumor_purity_l, function(d) {d["pur", "Median"]}) %>% unlist %>% append(NA)
    #     p2 = lapply(tumor_purity_l, function(d) {r = try(d["lambda.srt[1]", "Median"]); ifelse(class(r) == "try-error", NA, r)}) %>% unlist %>% append(NA)
    purity_dt = data.table(Chr="", Start="", End="", Ref="", Alt="", data.table(t(p1)), Gene.refGene=c("Tumor_Purity1"), Func.refGene="", ExonicFunc.refGene="", AAChange.refGene="")

    names(purity_dt) = names(m_sub)
    m_sub = rbind(m_sub, purity_dt)

    m_sub
}


# oncogene = fread("../data/allAnnotatedVariants.txt")$Gene %>% unique
symbol_list = c("APC", "TP53", "KRAS", "BRAF", "AMER1", "CTNNB1", "FBXW7", "Tumor_Purity1", "Tumor_Purity2") %>% unique
lapply(1:length(files), function(i) {
    read_annoted_vaf(i)
    }) -> tab_l

lapply(tab_l, function(tab) {
    measure_var_names = grep("FRAC", names(tab), value = T)
    tab_long = melt(tab
        , id.vars = c("Gene.refGene", "Start", "Ref", "Alt", "ExonicFunc.refGene")
        , measure.vars = measure_var_names
        , variable.name = "tumor_site"
        , value.name = "VAF"
    )
    tab_long
    }) -> tab_l2

tab_dt = rbindlist(tab_l2)
write_tsv(tab_dt, "../data/multiregion_purbyes_estimateion.tsv")

