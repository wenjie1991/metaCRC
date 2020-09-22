# input   ../data/01_tumor_purity.tsv
# output  Figure.S8 
library(ggplot2)
library(data.table)

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


## PurBayes tumor purity estimation
# tumor_purity = fread("../data/01_tumor_purity.tsv")
tumor_purity = fread("../data/01_tumor_purity_syn.tsv")
dim(tumor_purity)

tumor_purity_l = tumor_purity[, .(personID = c(personID, personID), median = c(median_p, median_m), p5 = c(p5_p, p5_m), p95 = c(p95_p, p95_m), sampleType = factor(rep(c("Primary", "Metastasis"), each = length(personID)), levels = c("Primary", "Metastasis")))]

tumor_purity_l[grepl("ID", personID), datasets := "Local Exome Seq"]
tumor_purity_l[grepl("AMC", personID), datasets := "Lim_2015"]

wilcox.test(tumor_purity_l[sampleType == "Primary", median] - tumor_purity_l[sampleType == "Metastasis", median])

set.seed(1)
#+ fig.height = 7, fig.width = 7, dev="pdf"
ggplot(data = tumor_purity_l, aes(x = sampleType, y=median, color = personID)) + 
    geom_boxplot(aes(x = sampleType, y = median), color ="black",  alpha=0) + 
    geom_point(size = 4, position=position_dodge(width=0.5), alpha=0.4) +
    geom_errorbar(
        aes(ymin = p5, ymax = p95),
        width = 0.1,
        position=position_dodge(width=0.5), alpha=0.4) +
    ylim(0, 1) + 
    ylab("Tumor purity") + 
    xlab("Sample Type") + theme0

tumor_purity[, .(mean = mean(median_p), sd = sd(median_p))]
tumor_purity[, .(mean = mean(median_m), sd = sd(median_m))]

## Facet tumor purity estimation
tumor_purity_facet = fread("../data/tumor_purity_facet.tsv") %>% na.omit
dim(tumor_purity)

tumor_purity_l = tumor_purity_facet[, .(personID = c(personID, personID), median = c(purity_p, purity_m), sampleType = factor(rep(c("Primary", "Metastasis"), each = length(personID)), levels = c("Primary", "Metastasis")))]
dim(tumor_purity_l)

tumor_purity_l[grepl("ID", personID), datasets := "Local Exome Seq"]
tumor_purity_l[grepl("AMC", personID), datasets := "Lim_2015"]

wilcox.test(tumor_purity_l[sampleType == "Primary", median] - tumor_purity_l[sampleType == "Metastasis", median])

#+ fig.height = 7, fig.width = 7, dev="pdf"
set.seed(1)
ggplot(data = tumor_purity_l, aes(x = sampleType, y=median, color = personID)) + 
    geom_boxplot(aes(x = sampleType, y = median), color ="black",  alpha=0) + 
    geom_point(size = 4, position=position_dodge(width=0.5), alpha=0.4) +
    ylim(0, 1) + 
    ylab("Tumor purity") + 
    xlab("Sample Type") + theme0

tumor_purity_facet[, .(mean = mean(purity_p), sd = sd(purity_p))]
tumor_purity_facet[, .(mean = mean(purity_m), sd = sd(purity_m))]

## Correlation between two tumor purity estimation
tumor_purity = fread("../data/01_tumor_purity_syn.tsv")

tumor_purity_l_pb = tumor_purity[, .(personID = c(personID, personID), median = c(median_p, median_m), sampleType = factor(rep(c("Primary", "Metastasis"), each = length(personID)), levels = c("Primary", "Metastasis")))]

tumor_purity_facet = fread("../data/tumor_purity_facet.tsv")

tumor_purity_l_facet = tumor_purity_facet[, .(personID = c(personID, personID), median = c(purity_p, purity_m), sampleType = factor(rep(c("Primary", "Metastasis"), each = length(personID)), levels = c("Primary", "Metastasis")))]

tumor_purity_merged = merge(tumor_purity_l_facet, tumor_purity_l_pb, by = c("personID", "sampleType"), suffixes = c("_Facet", "_PurBayes"))

#+ fig.height = 7, fig.width = 7, dev="pdf"
ggplot(data = tumor_purity_merged) + aes(x = median_Facet, y = median_PurBayes) +
    geom_point(size = 2) + 
    geom_smooth(method = "lm", formula = y ~ x) + theme0 +
    xlim(0, 1) +
    ylim(0, 1) +
    xlab("FACETS Tumor purity") + 
    ylab("PurBayes Tumor purity")


dim(tumor_purity_merged)
with(tumor_purity_merged, cor.test(x = median_Facet, y = median_PurBayes, method = "spearman"))
with(tumor_purity_merged, cor.test(x = median_Facet, y = median_PurBayes, method = "pearson"))
