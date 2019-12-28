# input         ../data/multiregion_biopsies/03_annotation/
# output        Figure2C, Figure2D
library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(readr)
library(plyr)
library(parallel)
library(PurBayes)
library(parallel)


# VAF files
files = dir("../data/multiregion_biopsies/03_annotation/", pattern="tsv", full=T)
person = str_match(basename(files), "\\w+")[, 1]

calc_tumor_purity = function(af_ref, af_alt) {
    ## 
    #     af_ref = af_ref_sub
    #     af_alt = af_alt_sub
    ##

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

plot_line_graph = function(tab, normal_name, primary_name, metastasis_name) {
    names(tab) %<>% sub("\\.FRAC", "", .)
    tab_long = melt(tab
        , id.vars = c("Gene.refGene", "Start", "Ref", "Alt", "ExonicFunc.refGene")
        , measure.vars = c(normal_name, primary_name, metastasis_name)
        , variable.name = "tumor_site"
        , value.name = "VAF"
        )

    tab_long[, mut_site := paste0(Gene.refGene, ":",  "c.", Ref,  Start, Alt)]
    tab_long$VAF %<>% as.numeric

    ggplot(tab_long) +
        aes(x=tumor_site, y=VAF) +
        geom_line(aes(color=mut_site, group=mut_site)) +
        geom_point(aes(shape=ExonicFunc.refGene, color=mut_site), size=I(3.5)) 
}

#+ fig.width=12, fig.height=7, dev='pdf'
oncogene = fread("../data/allAnnotatedVariants.txt")$Gene %>% unique
symbol_list = c("APC", "TP53", "KRAS", "BRAF", "AMER1", "CTNNB1", "FBXW7", "Tumor_Purity1", "Tumor_Purity2", oncogene) %>% unique
# tab1 = read_annoted_vaf(1)
# tab2 = read_annoted_vaf(2)
# tab3 = read_annoted_vaf(3)
# tab4 = read_annoted_vaf(4)
# tab5 = read_annoted_vaf(5)


normal_name = "CHET40_ND"
primary_name = c(
    "CHET40_1D"
    , "CHET40_3D"
    , "CHET40_4D" , "CHET40_5D"
    , "CHET40_6D"
    )
metastasis_name = c(
    "CHET40_LiMe2"
    , "CHET40_LiMe3"
    )
tab1 = tab1[Gene.refGene %in% symbol_list]
plot_line_graph(tab1, normal_name, primary_name, metastasis_name)

normal_name = "CHET54_ND"
primary_name = c(
    "CHET54_1D"
    , "CHET54_3D"
    , "CHET54_4D"
    , "CHET54_5D"
    , "CHET54_6D"
    )
metastasis_name = c(
    "CHET54_Lime_1D"
    , "CHET54_Lime_2D"
    , "CHET54_Lime_3D"
    , "CHET54_Lime_4D"
    , "CHET54_Lime_5D"
    , "CHET54_Lime_6D"
    )
tab2 = tab2[Gene.refGene %in% symbol_list]
plot_line_graph(tab2, normal_name, primary_name, metastasis_name)


normal_name = "CHET9_ND"
primary_name = c(
    "CHET9_1D"
    , "CHET9_5D"
    , "CHET9_6D"
    , "CHET9_7D"
    )
metastasis_name = c(
    "CHET9_LiMe"
    , "CHET9_LiMe2"
    , "CHET9_LiMe4"
    )
tab3 = tab3[Gene.refGene %in% symbol_list]
plot_line_graph(tab3, normal_name, primary_name, metastasis_name)


normal_name = "CHET_58_ND"
primary_name = c(
    "CHET_58_1D"
    , "CHET_58_2D"
    , "CHET_58_5D"
    )
metastasis_name = c(
    "CHET_58_LIME_1D"
    , "CHET_58_LIME_2D"
    , "CHET_58_LIME_3D"
    , "CHET_58_LIME_4D"
    )
tab4 = tab4[Gene.refGene %in% symbol_list]
plot_line_graph(tab4, normal_name, primary_name, metastasis_name)


normal_name = "Ecol_523_ND"
primary_name = c(
    "Ecol_523_1D"
    , "Ecol_523_2D"
    )
metastasis_name = c(
    "Ecol_523_LiME1D"
    , "Ecol_523_LiME2D"
    )
tab5 = tab5[Gene.refGene %in% symbol_list]
plot_line_graph(tab5, normal_name, primary_name, metastasis_name)

