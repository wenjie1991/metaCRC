#input      ../data/02_early_stage_CC_mut_rate.tsv
#           ../data/02_meta_stage_CC_mut_rate.tsv
#           ../data/02_meta_CC_mut_rate.tsv
#           ../data/02_Primary_mut_rate.tsv
#           ../data/02_Metastasis_mut_rate.tsv
# output    Figure1A 


library(magrittr)
library(data.table)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(plyr)


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

#' # Read data
# TCGA
# dat_early_tcga = fread("../data/02_early_stage_TCGA_mut_rate.tsv")
# dat_meta_tcga = fread("../data/02_advanced_stage_TCGA_mut_rate.tsv")
# CC
dat_early_cc = fread("../data/02_early_stage_CC_mut_rate.tsv")
(n_early_cc = dat_early_cc[1, freq_n / freq * 100 %>% round])
dat_meta_cc = fread("../data/02_meta_stage_CC_mut_rate.tsv")
(n_meta_cc = dat_meta_cc[1, freq_n / freq * 100 %>% round])
dat_liver_cc = fread("../data/02_meta_CC_mut_rate.tsv")
(n_liver_cc = dat_liver_cc[1, freq_n / freq * 100 %>% round])
# Trios
dat_trios_p = fread("../data/02_Primary_mut_rate.tsv")
dat_trios_m = fread("../data/02_Metastasis_mut_rate.tsv")
(n_trios = dat_trios_p[1, freq_n / freq * 100 %>% round])


#' # Combind data 
dat_combind = merge(dat_early_cc, dat_meta_cc, all=T, by="symbol", suffixes=c("_early_cc", "_meta_cc")) %>% 
    merge(dat_trios_m, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_trios_m", "freq"="freq_trios_m")) %>%
    merge(dat_trios_p, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_trios_p", "freq"="freq_trios_p")) %>% 
    merge(dat_liver_cc, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_liver_cc", "freq"="freq_liver_cc"))
    #     merge(dat_early_tcga, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_early_tcga", "freq"="freq_early_tcga")) %>% 
    #     merge(dat_meta_tcga, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_meta_tcga", "freq"="freq_meta_tcga"))
dat_combind[dat_combind %>% is.na] = 0

#' # Output Table
dat_combind[freq_trios_m > 3 | freq_trios_p > 3, .(symbol, freq_early_cc, freq_meta_cc, freq_liver_cc, freq_trios_p, freq_trios_m)]
write_tsv(dat_combind, "../data/03_mutation_frequency.tsv")

#' # Statistics test
dat_combind_sub = dat_combind[freq_trios_m > 3 | freq_trios_p > 3]
fisher_result = sapply(1:nrow(dat_combind_sub), function(i) {
    m = dat_combind_sub[i, .(
        freq_n_early_cc, n_early_cc - freq_n_early_cc, 
        freq_n_meta_cc, n_meta_cc - freq_n_meta_cc,
        freq_n_liver_cc, n_liver_cc - freq_n_liver_cc,
        freq_n_trios_p, n_trios - freq_n_trios_p,
        freq_n_trios_m, n_trios - freq_n_trios_m)] %>% unlist %>% matrix(nrow=2)
    c(
        all = fisher.test(m)$p.value, 
        tP_tM = fisher.test(m[, c(4, 5)])$p.value,
        cP_cM = fisher.test(m[, c(2, 3)])$p.value,
        cE_tP = fisher.test(m[, c(1, 4)])$p.value,
        cE_tM = fisher.test(m[, c(1, 5)])$p.value,
        cE_cP = fisher.test(m[, c(1, 2)])$p.value,
        cE_cM = fisher.test(m[, c(1, 3)])$p.value
        )
    })
fisher_result = data.frame(t(fisher_result))

dat_combind_sub = cbind(dat_combind_sub, fisher_result)

mut_Freq_test_result = dat_combind_sub[,.(symbol, freq_n_early_cc, freq_early_cc, freq_meta_cc, freq_liver_cc, freq_trios_p, freq_trios_m)] %>% cbind(fisher_result)
write_tsv(mut_Freq_test_result, "../data/03_mut_Freq_test.tsv")


#' # Draw graph
#+ mutation Freq by stage,  fig.height=7, fig.width=15, dev='pdf'
# dat_combind_long = dat_combind[freq_trios_m > 5 | freq_trios_p > 5 | freq_liver_cc > 5 | freq_early_cc > 5 | freq_meta_cc > 5, .(symbol, freq_early_cc, freq_meta_cc, freq_liver_cc, freq_trios_p, freq_trios_m)] %>% melt
dat_combind_long = dat_combind[freq_trios_m > 5 | freq_trios_p > 5, .(symbol, freq_early_cc, freq_meta_cc, freq_liver_cc, freq_trios_p, freq_trios_m)] %>% melt
dat_combind_long %<>% rename(c("variable" = "group", "value" = "mutation_freq"))
dat_combind_long$symbol %<>% factor(levels = dat_combind_long[, .(s = sum(mutation_freq)), by=symbol][order(s, decreasing=T), symbol])

colors = brewer.pal(5, "Blues")[3:5]
colors[4:5] =  brewer.pal(3, "Greens")[c(2, 3)]
ggplot(data = dat_combind_long) +
    aes(x=symbol, y=mutation_freq, label=paste0(round(mutation_freq, 1), "%"), fill=group) +
    coord_cartesian(ylim = c(0, max(dat_combind_long$mutation_freq) * 1.1)) +
    geom_bar(data = dat_combind_long, position = "dodge", stat = "identity") +
    geom_text(position=position_dodge(width=.9), vjust=0.7, hjust=-0.2, angle=90, size=4.5) +
    scale_fill_manual(values=colors) + theme
