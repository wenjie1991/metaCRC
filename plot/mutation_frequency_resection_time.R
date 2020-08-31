# input     ../data/02_Primary_mut_rate_syn.tsv
#           ../data/02_Primary_mut_rate_sub.tsv
#           ../data/02_Metastasis_mut_rate_syn.tsv
#           ../data/02_Metastasis_mut_rate_sub.tsv
# output    FigureS5

library(magrittr)
library(data.table)
library(readr)
library(ggplot2)
library(RColorBrewer)

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
# Trios
dat_trios_p_syn = fread("../data/02_Primary_mut_rate_syn.tsv")
dat_trios_p_sub = fread("../data/02_Primary_mut_rate_sub.tsv")
dat_trios_m_syn = fread("../data/02_Metastasis_mut_rate_syn.tsv")
dat_trios_m_sub = fread("../data/02_Metastasis_mut_rate_sub.tsv")

(n_syn_trios = dat_trios_p_syn[1, freq_n / freq * 100])
(n_sub_trios = dat_trios_p_sub[1, freq_n / freq * 100])

#' # Combind data 
dat_combind = merge(dat_trios_p_syn, dat_trios_p_sub, all=T, by="symbol", suffixes=c("_p_syn", "_p_sub")) %>% 
    merge(dat_trios_m_syn, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_m_syn", "freq"="freq_m_syn")) %>%
    merge(dat_trios_m_sub, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_m_sub", "freq"="freq_m_sub"))
    #     merge(dat_early_tcga, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_early_tcga", "freq"="freq_early_tcga")) %>% 
    #     merge(dat_meta_tcga, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_meta_tcga", "freq"="freq_meta_tcga"))
dat_combind[dat_combind %>% is.na] = 0

#' # Output Table
dat_combind[freq_p_syn > 5 | freq_p_sub > 5 | freq_m_syn > 5| freq_m_sub > 5, .(symbol, freq_p_syn, freq_p_sub, freq_m_syn, freq_m_sub)]
write_tsv(dat_combind, "../data/03_mutation_frequency_resection_time.tsv")

#' # Statistics test
dat_combind_sub = dat_combind[freq_p_syn > 10 | freq_p_sub > 10 | freq_m_syn > 10 | freq_m_sub > 10]
fisher_result = sapply(1:nrow(dat_combind_sub), function(i) {
    m = dat_combind_sub[i, .(
        freq_n_p_syn, n_syn_trios - freq_n_p_syn, 
        freq_n_p_sub, n_sub_trios - freq_n_p_sub,
        freq_n_m_syn, n_syn_trios - freq_n_m_syn,
        freq_n_m_sub, n_sub_trios - freq_n_m_sub)] %>% unlist %>% matrix(nrow=2)
    c(
        all = fisher.test(m)$p.value, 
        pSyn_pSub = fisher.test(m[, c(1, 2)])$p.value,
        mSyn_mSub = fisher.test(m[, c(3, 4)])$p.value,
        pSyn_mSyn = fisher.test(m[, c(1, 3)])$p.value,
        pSub_mSub = fisher.test(m[, c(2, 4)])$p.value
        )
    })

fisher_result = data.frame(t(fisher_result))
dat_combind_sub = cbind(dat_combind_sub, fisher_result)

mut_Freq_test_result = dat_combind_sub
write_tsv(mut_Freq_test_result, "../data/03_mut_Freq_test_resection_time.tsv")


#' # Draw graph
#+ mutation Freq by resection time,  fig.height=7, fig.width=15, dev='pdf'
# dat_combind_long = dat_combind[freq_trios_m > 5 | freq_trios_p > 5 | freq_liver_cc > 5 | freq_early_cc > 5 | freq_meta_cc > 5, .(symbol, freq_early_cc, freq_meta_cc, freq_liver_cc, freq_trios_p, freq_trios_m)] %>% melt
dat_combind_long = dat_combind_sub[, .(symbol, freq_p_syn, freq_p_sub, freq_m_syn, freq_m_sub)] %>% melt
dat_combind_long %<>% rename(c("variable" = "group", "value" = "mutation_freq"))
dat_combind_long$symbol %<>% factor(levels = dat_combind_long[, .(s = sum(mutation_freq)), by=symbol][order(s, decreasing=T), symbol])

colors = brewer.pal(5, "Blues")[3:5]
colors[3:4] =  brewer.pal(4, "Greens")[c(2, 4)]
ggplot(data = dat_combind_long) +
    aes(x=symbol, y=mutation_freq, label=paste0(round(mutation_freq, 1), "%"), fill=group) +
    coord_cartesian(ylim = c(0, max(dat_combind_long$mutation_freq) * 1.1)) +
    geom_bar(data = dat_combind_long, position = "dodge", stat = "identity") +
    geom_text(position=position_dodge(width=.9), vjust=0.7, hjust=-0.2, angle=90, size=4) +
    scale_fill_manual(values=colors)  + theme
