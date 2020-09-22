# input         ../data/result_depth50/result_dt.tsv
#               ../data/result_depth50/result_dt.tsv
#               ../data/phenotype.csv
#               ../data/00_mean_coverage.tsv
#               ../data/VAP_real_data/patient_metrics.csv
# output        Figure5A, Figure5B, Figure5C

library(ggplot2)
library(magrittr)
library(data.table)


library(RColorBrewer)
colors = brewer.pal(5, "Blues")[c(3, 5)]
colors[3:5] =  brewer.pal(5, "Greens")[c(3:5)]

dat = fread("../data/result_depth50/result_dt.tsv")
dat[, SampleType := revalue(subA, c("maf1" = "primary", "maf2" = "primary", "maf3" = "primary", "maf4" = "primary", "maf5" = "metastasis", "maf6" = "metastasis", "maf7" = "metastasis", "maf8" = "metastasis")) %>% paste("-", revalue(subB, c("maf1" = "primary", "maf2" = "primary", "maf3" = "primary", "maf4" = "primary", "maf5" = "metastasis", "maf6" = "metastasis", "maf7" = "metastasis", "maf8" = "metastasis")), sep="")]
dat[, KSD := as.numeric(KSD)]
dat[, FST := as.numeric(FST)]

dat_line = dat[, .(
    FST_mean = mean(FST)
    , FST_sd = sd(FST) / sqrt(.N) * 1.96
    , KSD_mean = mean(KSD)
    , KSD_sd = sd(KSD) / sqrt(.N) * 1.96
    ), by = .(selection, SampleType, meta_init, meta_cell_n)]


## Line plot
#+ KSD line, fig.height=6, fig.width=9, dev="pdf"
ggplot(dat_line[selection == 0 & SampleType == "primary-metastasis"]) + aes(x = factor(meta_init), y = KSD_mean, color = factor(meta_cell_n)) + geom_pointrange(aes(ymin = KSD_mean - KSD_sd, ymax = KSD_mean + KSD_sd), alpha=0.7, size=1.1) + geom_line(aes(group=factor(meta_cell_n))) + scale_color_manual(values=colors) + labs(y = "KSD", x = "Primary tumor cell number", color = "Metastasis\nCell Number")

#+ FST line, fig.height=6, fig.width=9, dev="pdf"
ggplot(dat_line[selection == 0 & SampleType == "primary-metastasis"]) + aes(x = factor(meta_init), y = FST_mean, color = factor(meta_cell_n)) + geom_pointrange(aes(ymin = FST_mean - FST_sd, ymax = FST_mean + FST_sd), alpha=0.7, size=1.1) + geom_line(aes(group=factor(meta_cell_n))) + scale_color_manual(values=colors) + labs(y = "FST", x = "Primary tumor cell number", color = "Metastasis\nCell Number")

#+ KSD density, fig.height=6, fig.width=10, dev="pdf"
ggplot(dat[SampleType == "primary-metastasis" & selection == 0 & meta_cell_n %in% c(1, 10, 100, 1000)], aes(x = KSD)) + 
    geom_density(aes(fill = factor(meta_cell_n)), alpha=0.5) + scale_fill_manual(values=colors[c(1, 2, 3, 5)]) +
    theme(axis.text.x = element_text( margin = margin(t = 2, unit = "mm") , angle = 0, vjust = 1, size = 25, hjust = 0.5)) +
    labs(fill = "Metastasis\nCell Number", y = "Density")

#+ FST density, fig.height=6, fig.width=10, dev="pdf"
ggplot(dat[SampleType == "primary-metastasis" & selection == 0 & meta_cell_n %in% c(1, 10, 100, 1000)], aes(x = FST)) + 
    geom_density(aes(fill = factor(meta_cell_n)), alpha=0.5) + scale_fill_manual(values=colors[c(1, 2, 3, 5)]) +
    theme(axis.text.x = element_text( margin = margin(t = 2, unit = "mm") , angle = 0, vjust = 1, size = 25, hjust = 0.5)) +
    labs(fill = "Metastasis\nCell Number", y = "Density")
## Real data
library(stringr)


cd = fread("../data/phenotype.csv")
cd[grepl("^20", pathoID), patientID := paste0("ID", substr(pathoID, 3, 9))]


coverage = fread("../data/00_mean_coverage.tsv")
unqualify_sample = coverage[mean < 50 & !grepl("sn$", sample_name), sample_name] %>% str_match(., "(^\\w+)_") %>% extract(, 2) %>% unique

real_dat = fread("../data/VAP_real_data/patient_metrics.csv")[!(patientID %in% c("ID0825371", "ID1119999"))]
dat_sub1 = real_dat[grepl("^(AMC)|(ID11)", patientID)][, subset := paste0("1.All(", .N, ")")][!is.na(KSD)]
dat_sub2 = dat_sub1[!(grepl("^ID", patientID) & patientID %in% cd[Both_treated != "Chemonaive", patientID])][, subset := paste0("2.Remove treated(", .N, ")")]
dat_sub3 = dat_sub2[!patientID %in% unqualify_sample][, subset := paste0("3.Remove low depth(", .N, ")")]

dat_for_plot = rbind(dat_sub1, dat_sub2, dat_sub3)

#+ real KSD QC, fig.height=6, fig.width=10, dev="pdf"
ggplot(dat_for_plot) + aes(x = KSD) + geom_density(aes(fill = subset), alpha = 0.5) + 
    theme(axis.text.x = element_text( margin = margin(t = 2, unit = "mm") , angle = 0, vjust = 1, size = 25, hjust = 0.5)) +
    labs(fill = "QC process", y = "Density")



#+ combine KSD, fig.height=6, fig.width=10, dev="pdf"
dat_for_plot = rbind(dat[SampleType == "primary-metastasis" & selection == 0 & meta_cell_n %in% c(1, 10, 100, 1000), .(KSD, meta_cell_n)],
    dat_sub3[, .(KSD, meta_cell_n = "Real KSD")])

ggplot(dat_for_plot, aes(x = KSD)) + 
    geom_density(aes(fill = factor(meta_cell_n)), alpha=0.5) + scale_fill_manual(values=c(colors[c(1, 2, 3, 5)], "red")) +
    theme(axis.text.x = element_text( margin = margin(t = 2, unit = "mm") , angle = 0, vjust = 1, size = 25, hjust = 0.5)) +
    labs(fill = "Metastasis\nCell Number", y = "Density")

ggplot(dat_for_plot, aes(x = meta_cell_n, y = KSD, fill = meta_cell_n)) + 
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values=c(colors[c(1, 2, 3, 5)], "red")) +
    theme(axis.text.x = element_text( margin = margin(t = 2, unit = "mm") , angle = 0, vjust = 1, size = 25, hjust = 0.5)) +
    labs(fill = "Metastasis\nCell Number", y = "KSD", x = "Metastasis Cell Number")


wilcox.test(dat_for_plot[meta_cell_n == "1", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])
wilcox.test(dat_for_plot[meta_cell_n == "10", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])
wilcox.test(dat_for_plot[meta_cell_n == "100", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])
wilcox.test(dat_for_plot[meta_cell_n == "1000", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])

wilcox.test(dat_for_plot[meta_cell_n == "1", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])
wilcox.test(dat_for_plot[meta_cell_n == "10", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])
wilcox.test(dat_for_plot[meta_cell_n == "100", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])
wilcox.test(dat_for_plot[meta_cell_n == "1000", KSD], dat_for_plot[meta_cell_n == "Real KSD", KSD])
