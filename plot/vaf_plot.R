# input         ../data/01_trios_somaticMutation_raw.txt
#               ../data/01_trios_somaticMutation.txt")
#               ../data/tumor_purity_facet.tsv
# output        Figure2.A, Figure2.B
#               ../data/02_mutDT.tsv

library(magrittr)
library(data.table)
library(ggplot2)
library(plyr)
library(readr)

tumorPurityThreshold = 0.10
normalPurityThreshold = 0.01

theme <- theme_bw() + theme(
    text = element_text(size = 20),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 60, vjust = 1, size = 15, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "top",
    ) 

plot_AF_line = function(geneSymbol) {
    geneSymbol = "TP53"
    data = deltaAFlong_raw[Gene.refGene %in% geneSymbol & loci %in% deltaAFlong$loci, .(sampleType = AF_Type, AF, sampleMut = paste0(personID, loci), source)]
    p = ggplot(data=data, aes(x = sampleType %>% revalue(c(AF.N = "Normal", AF.M = "Metastasis", AF.C = "Primary")), y = AF)) +
        geom_line(aes(group = sampleMut, color = source)) +
        geom_point(aes(color = source)) +
        theme + theme(legend.position = "right", axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5)) + 
        labs(title = paste0(geneSymbol, " mutation AF")) + xlab("Tumor Site") + ylab("VAF") +
        ylim(0, 1)
    p
}

#+ Read in data
annTabs_raw <- fread("../data/01_trios_somaticMutation_raw.txt")
annTabs2.0 <- fread("../data/01_trios_somaticMutation.txt")
msk_impact_341 = readLines("../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv")
annTabs2.1 <- annTabs2.0[(metastasis_site %in% c("Liver")) | grepl("AMC", personID)]
annTabs2.2 <- annTabs2.1[MS == "MSS"]
annTabs2 = annTabs2.2[Gene.refGene %in% msk_impact_341]
n_metastasis = annTabs_raw$personID %>% table %>% length
n_right_metastasis = annTabs_raw[primary_site == "right"]$personID %>% table %>% length
n_left_metastasis = annTabs_raw[primary_site == "left"]$personID %>% table %>% length
n_syn_metastasis = annTabs_raw[Resection_timing == "Concurrent"]$personID %>% table %>% length
n_sub_metastasis = annTabs_raw[Resection_timing == "Subsequent"]$personID %>% table %>% length
n_naive_metastasis = annTabs_raw[Both_treated == "Chemonaive"]$personID %>% table %>% length
n_both_metastasis = annTabs_raw[Both_treated == "Both_treated"]$personID %>% table %>% length
n_meta_metastasis = annTabs_raw[Both_treated == "Metastasis_treated"]$personID %>% table %>% length

# Cancer cell frequency
ccf_d = fread("../data/tumor_purity_facet.tsv") %>% na.omit

#' # Count the sample size
(annTabs_raw[, .(personID, primary_site)] %>% unique)[, primary_site] %>% table(useNA = "ifany")
(annTabs_raw[, .(personID, primary_site)] %>% unique)[, primary_site] %>% table(useNA = "ifany") %>% sum
(annTabs_raw[, .(personID, Both_treated)] %>% unique)[, Both_treated] %>% table(useNA = "ifany")
(annTabs_raw[, .(personID, Resection_timing)] %>% unique)[, Resection_timing] %>% table(useNA = "ifany")
annTabs_raw[, unique(personID) %>% sub("\\d+", "", .) %>% table]
annTabs_raw[, unique(personID)]


somaticMutation <- annTabs2[(C.A / (C.A + C.R + 1) > tumorPurityThreshold & C.A > 2) | (M.A / (M.A + M.R + 1) > tumorPurityThreshold & M.A > 2) 
    , .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]

#' ## Primary Mutation
somaticMutationP <- annTabs2[(C.A / (C.A + C.R + 1) > tumorPurityThreshold & C.A > 2), .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]


somaticMutationM <- annTabs2[(M.A / (M.A + M.R + 1) > tumorPurityThreshold & M.A > 2), .(Gene.refGene, personID)][, .(n=length(unique(personID))), by=Gene.refGene][order(n, decreasing=T)]

somaticMutationP[, group := "primary"]
somaticMutationM[, group := "metastasis"]
somaticMutationPM <- rbind(somaticMutationP, somaticMutationM)
somaticMutationPM$Gene.refGene %<>% factor(
    levels = somaticMutationPM[, .(s=sum(n)), by=Gene.refGene][order(s, decreasing=T), Gene.refGene])
somaticMutation$Gene.refGene %<>% factor(levels=somaticMutation$Gene.refGene)

#' ### Allelle Frequency 
deltaAF_all <- annTabs2[, .(
    personID
    , loci = paste0(CHROM, ":", POS, REF, "/", ALT)
    , cosmic70 = cosmic70
    , AAChange.refGene = AAChange.refGene
    , source= personID %>% sub("AMC\\d+", "AMC1", .) %>%  gsub("\\d", "+", .)
    , AF.N = N.A / (N.A + N.R + 1) %>% as.numeric
    , AF.M = M.A / (M.A + M.R + 1) %>% as.numeric
    , AF.C = C.A / (C.A + C.R + 1) %>% as.numeric
    , depth.N = N.A + N.R %>% as.numeric
    , depth.M = M.A + M.R %>% as.numeric
    , depth.C = C.A + C.R %>% as.numeric
    , Gene.refGene)] %>% data.table
AF_all = deltaAF_all
AF_all= merge(AF_all, ccf_d, by = "personID", all.x = T)
AF_all[is.na(purity_p), purity_p := 1]
AF_all[is.na(purity_m), purity_m := 1]
AF_all[, deltaAF := (AF.M) - (AF.C)]
AF_all[, deltaAF_corrected := (AF.M / purity_m) - (AF.C / purity_p)]

deltaAF = AF_all
deltaAF$source %<>% revalue(c("ID+++++++"="Local", "AMC+"="Lim_2015", "EV-+++"="Brannon_2014", "P-+++++++"="Yaeger_2018"))

deltaAFlong_raw <- melt(deltaAF, id.vars=c("personID", "deltaAF",  "deltaAF_corrected", "Gene.refGene", "source", "loci")
    , measure.vars=c("AF.M", "AF.C", "AF.N")
    , value.name="AF", variable="AF_Type")
deltaAFlong_raw$AF_Type %<>% factor(levels = c("AF.N", "AF.C", "AF.M"))


selected_genes = c("APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7", "AMER1", "BRAF", "ATM")
deltaAFlong = deltaAFlong_raw[AF_Type == "AF.C", .(AF = max(AF), deltaAF = deltaAF[which.max(AF)], deltaAF_corrected = deltaAF_corrected[which.max(AF)], loci=loci[which.max(AF)]), by=c("Gene.refGene", "source", "personID", "AF_Type")]
deltaAFlong = deltaAFlong[Gene.refGene %in% selected_genes]
deltaAFlong$Gene.refGene %<>% factor(levels = selected_genes)

#' ## VAF difference between Primary and metastasis 
#+ Delta Allele Frequency,  fig.height=7, fig.width=13, dev='svg'
deltaAFlong_draw = deltaAFlong[, .(Gene.refGene, deltaAF, source)] %>% unique
deltaAFlong_draw$Gene.refGene %<>% factor(levels = deltaAFlong_draw[, .N, by=Gene.refGene][order(N, decreasing=T), Gene.refGene])
ggplot(data= deltaAFlong_draw %>% unique, aes(x=Gene.refGene, deltaAF)) + 
    geom_boxplot(alpha=0.2) + 
    geom_jitter(width=0.3, aes(color=source), alpha=0.3, size=3) + 
    xlab("Gene") + ylab("VAF_Metastasis - VAF_Primary") + theme

#+ statistics Test
deltaAFlong_draw[, .(p = wilcox.test(deltaAF)$p.value, n = .N, median = median(deltaAF)), by = Gene.refGene]


#' ## VAF difference between Primary and metastasis with correction
#+ Delta Allele Frequency with correction,  fig.height=7, fig.width=13, dev='svg'
deltaAFlong_draw = deltaAFlong[, .(Gene.refGene, deltaAF_corrected, deltaAF, source)] %>% unique
deltaAFlong_draw$Gene.refGene %<>% factor(levels = deltaAFlong_draw[, .N, by=Gene.refGene][order(N, decreasing=T), Gene.refGene])
ggplot(data= deltaAFlong_draw %>% unique, aes(x=Gene.refGene, deltaAF_corrected)) + 
    geom_boxplot(alpha=0.2) + 
    geom_jitter(width=0.3, aes(color=source), alpha=0.3, size=3) + 
    xlab("Gene") + ylab("VAF_Metastasis - VAF_Primary") + theme

deltaAFlong_draw[, .(p = wilcox.test(deltaAF_corrected)$p.value, n = .N, median = median(deltaAF_corrected)), by = Gene.refGene]
deltaAFlong_draw[, .(p = wilcox.test(deltaAF)$p.value, n = .N, median = median(deltaAF)), by = Gene.refGene]

deltaAFlong_draw[!grepl("Lim|Local",source) , .(p = wilcox.test(deltaAF_corrected)$p.value, n = .N, median = median(deltaAF_corrected)), by = Gene.refGene]
deltaAFlong_draw[!grepl("Lim|Local",source), .(p = wilcox.test(deltaAF)$p.value, n = .N, median = median(deltaAF)), by = Gene.refGene]

#' ### Allele Frequency of small Wilcox P value
#+ AF by sample type,  fig.height=7, fig.width=10, dev='svg'
plot_AF_line("TP53")

#' ### Allele Frequency v.s. depth for TP53
#+ AF v.s. depth,  fig.height=7, fig.width=10, dev='svg'
d = deltaAF[Gene.refGene == "TP53" & loci %in% deltaAFlong$loci, .(source, AF.M, AF.C, depth.M, depth.C)]
d_long = d[, .(
    source = c(source, source)
    , AF = c(AF.M, AF.C)
    , depth = c(depth.M, depth.C)
    , sample_type = rep(c("Metastasis", "Primary"), each = length(source))
    )]
n_primary = d_long[sample_type == "Primary"] %>% nrow
n_metastasis = d_long[sample_type == "Metastasis"] %>% nrow

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
library(RColorBrewer)
colors = c(brewer.pal(5, "Blues")[5], brewer.pal(4, "Greens")[4])
g = ggplot(d_long) + aes(x = depth, y = AF, color = sample_type)
p = g + geom_point(alpha = 0.7) + theme0 
library(stringr)
color_label = c(
    str_glue("Primary Tumor (n={n_primary})"), 
    str_glue("Metastasis Tumor (n={n_metastasis})")
    )
p + xlab("Depth") + ylab("VAF") + labs(color = "Sample Type") + scale_color_manual(labels = color_label, values = colors)

g = ggplot(d_long) + aes(x = sample_type, y = depth) 
g + geom_violin() + theme0 + labs(x = "Sample Type", y = "Depth")

with(deltaAF_all[Gene.refGene == "TP53"], cor.test(AF.M, depth.M, method = "spearman"))
dim(deltaAF_all[Gene.refGene == "TP53"])
with(deltaAF_all[Gene.refGene == "TP53"], cor.test(AF.C, depth.C, method = "spearman"))
with(deltaAF_all[Gene.refGene == "TP53"], cor.test(c(AF.M,AF.C), c(depth.M, depth.C), method = "spearman"))

(d_long[sample_type == "Primary", depth] - d_long[sample_type == "Metastasis", depth]) %>% median
wilcox.test(d_long[sample_type == "Primary", depth] - d_long[sample_type == "Metastasis", depth])
