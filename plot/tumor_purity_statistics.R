# input   ../data/01_tumor_purity.tsv
# output  Figure.S8 
library(ggplot2)
library(data.table)

tumor_purity = fread("../data/01_tumor_purity.tsv")[! (personID %in% c("ID0825371", "ID119999"))]

tumor_purity_l = tumor_purity[, .(personID = c(personID, personID), median = c(median_p, median_m), p5 = c(p5_p, p5_m), p95 = c(p95_p, p95_m), sampleType = factor(rep(c("Primary", "Metastasis"), each = length(personID)), levels = c("Primary", "Metastasis")))]

tumor_purity_l[grepl("ID", personID), datasets := "Local Exome Seq"]
tumor_purity_l[grepl("AMC", personID), datasets := "Lim_2015"]

wilcox.test(tumor_purity_l[sampleType == "Primary", median] - tumor_purity_l[sampleType == "Metastasis", median])

set.seed(1)
ggplot(data = tumor_purity_l, aes(x = sampleType, y=median, color = personID)) + 
    geom_boxplot(aes(x = sampleType, y = median), color ="black",  alpha=0) + 
    geom_point(size = 4, position=position_dodge(width=0.5), alpha=0.4) +
    geom_errorbar(
        aes(ymin = p5, ymax = p95),
        width = 0.1,
        position=position_dodge(width=0.5), alpha=0.4) +
    ylab("Tumor purity") + 
    xlab("Sample Type") +
