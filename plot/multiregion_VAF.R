# input         ../data/multiregion_biopsies/03_annotation/
#               ../data/multiregion_biopsies_facet_purity.tsv
#               ../data/multiregion_purbyes_estimateion.tsv
# output        Figure2C, Figure2D

library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(readr)
library(plyr)


# VAF files
files = dir("../data/multiregion_biopsies/03_annotation/", pattern="tsv", full=T)
person = str_match(basename(files), "\\w+")[, 1]

plot_line_graph = function(tab_long, normal_name, primary_name, metastasis_name) {
    
   
    tab_long[, mut_site := paste0(Gene.refGene, ":",  "c.", Ref,  Start, Alt)]
    tab_long$VAF %<>% as.numeric
    tab_long$ExonicFunc.refGene %<>% (function(x) { x[is.na(x)] = ""; x })
    tab_long$tumor_site %<>% factor(levels = c(normal_name, primary_name, metastasis_name))

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

    ggplot(tab_long) +
        aes(x=tumor_site, y=VAF) +
        geom_line(aes(color=mut_site, group=mut_site)) +
        geom_point(aes(shape=ExonicFunc.refGene, color=mut_site), size=I(3.5)) + theme0
}

# oncogene = fread("../data/allAnnotatedVariants.txt")$Gene %>% unique
symbol_list = c("APC", "TP53", "KRAS", "BRAF", "AMER1", "CTNNB1", "FBXW7", "Tumor_Purity1", "Tumor_Purity2") %>% unique

tab_dt = fread("../data/multiregion_purbyes_estimateion.tsv")
tab_dt$tumor_site %<>% sub("\\.FRAC", "", .)
purity2 = fread("../data/multiregion_biopsies_facet_purity.tsv")
purity2$Gene.refGene = "Tumor_Purity2"
names(purity2)[1:2] = c("tumor_site", "VAF")

tab_dt = merge(tab_dt, purity2, all = T)

#+ fig.width=12, fig.height=7, dev='pdf'


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
samples = c(normal_name, primary_name, metastasis_name)
tab = tab_dt[Gene.refGene %in% symbol_list & tumor_site %in% samples]
plot_line_graph(tab, normal_name, primary_name, metastasis_name)


normal_name = "Ecol_523_ND"
primary_name = c(
    "Ecol_523_1D"
    , "Ecol_523_2D"
    )
metastasis_name = c(
    "Ecol_523_LiME1D"
    , "Ecol_523_LiME2D"
    )
samples = c(normal_name, primary_name, metastasis_name)
tab = tab_dt[Gene.refGene %in% symbol_list & tumor_site %in% samples]
plot_line_graph(tab, normal_name, primary_name, metastasis_name)

# normal_name = "CHET40_ND"
# primary_name = c(
#     "CHET40_1D"
#     , "CHET40_3D"
#     , "CHET40_4D" , "CHET40_5D"
#     , "CHET40_6D"
#     )
# metastasis_name = c(
#     "CHET40_LiMe2"
#     , "CHET40_LiMe3"
#     )
# tab1 = tab1[Gene.refGene %in% symbol_list]
# plot_line_graph(tab1, normal_name, primary_name, metastasis_name)


# normal_name = "CHET9_ND"
# primary_name = c(
#     "CHET9_1D"
#     , "CHET9_5D"
#     , "CHET9_6D"
#     , "CHET9_7D"
#     )
# metastasis_name = c(
#     "CHET9_LiMe"
#     , "CHET9_LiMe2"
#     , "CHET9_LiMe4"
#     )
# tab3 = tab3[Gene.refGene %in% symbol_list]
# plot_line_graph(tab3, normal_name, primary_name, metastasis_name)
# 
# 
# normal_name = "CHET_58_ND"
# primary_name = c(
#     "CHET_58_1D"
#     , "CHET_58_2D"
#     , "CHET_58_5D"
#     )
# metastasis_name = c(
#     "CHET_58_LIME_1D"
#     , "CHET_58_LIME_2D"
#     , "CHET_58_LIME_3D"
#     , "CHET_58_LIME_4D"
#     )
# tab4 = tab4[Gene.refGene %in% symbol_list]
# plot_line_graph(tab4, normal_name, primary_name, metastasis_name)
