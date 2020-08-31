# input     ../data/02_Primary_mut_rate_both.tsv
#           ../data/02_Primary_mut_rate_naive.tsv
#           ../data/02_Primary_mut_rate_meta.tsv
#           ../data/02_Metastasis_mut_rate_both.tsv
#           ../data/02_Metastasis_mut_rate_naive.tsv
#           ../data/02_Metastasis_mut_rate_meta.tsv
# output    FigureS4

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
dat_trios_p_both = fread("../data/02_Primary_mut_rate_both.tsv")
dat_trios_p_naive = fread("../data/02_Primary_mut_rate_naive.tsv")
dat_trios_p_meta = fread("../data/02_Primary_mut_rate_meta.tsv")
dat_trios_m_both = fread("../data/02_Metastasis_mut_rate_both.tsv")
dat_trios_m_naive = fread("../data/02_Metastasis_mut_rate_naive.tsv")
dat_trios_m_meta = fread("../data/02_Metastasis_mut_rate_meta.tsv")

(n_both_trios = dat_trios_p_both[1, freq_n / freq * 100 %>% round])
(n_naive_trios = dat_trios_p_naive[1, freq_n / freq * 100 %>% round])
(n_meta_trios = dat_trios_p_meta[1, freq_n / freq * 100 %>% round])

#' # Combind data 
dat_combind = merge(dat_trios_p_both, dat_trios_p_naive, all=T, by="symbol", suffixes=c("_p_both", "_p_naive")) %>% 
    merge(dat_trios_p_meta, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_p_meta", "freq"="freq_p_meta")) %>%
    merge(dat_trios_m_both, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_m_both", "freq"="freq_m_both")) %>% 
    merge(dat_trios_m_naive, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_m_naive", "freq"="freq_m_naive")) %>% 
    merge(dat_trios_m_meta, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_m_meta", "freq"="freq_m_meta"))
dat_combind[dat_combind %>% is.na] = 0

#' # Output Table
dat_combind[freq_p_both > 10 | freq_p_naive > 10 | freq_p_meta > 10 | freq_m_both > 10| freq_m_naive > 10 | freq_m_meta > 10, .(symbol, freq_p_both, freq_p_naive, freq_p_meta, freq_m_both, freq_m_naive, freq_m_meta)]
write_tsv(dat_combind, "../data/03_mutation_frequency_therapy.tsv")

#' # Statistics test
dat_combind_sub = dat_combind[freq_p_both > 10 | freq_p_naive > 10 | freq_p_meta > 10 | freq_m_both > 10| freq_m_naive > 10 | freq_m_meta > 10 ]
fisher_result = sapply(1:nrow(dat_combind_sub), function(i) {
    m = dat_combind_sub[i, .(
        freq_n_p_naive, n_naive_trios - freq_n_p_naive,
        freq_n_p_both, n_both_trios - freq_n_p_both, 
        freq_n_p_meta, n_meta_trios - freq_n_p_meta,
        freq_n_m_naive, n_naive_trios - freq_n_m_naive,
        freq_n_m_both, n_both_trios - freq_n_m_both,
        freq_n_m_meta, n_meta_trios - freq_n_m_meta)] %>% unlist %>% matrix(nrow=2)
    c(
        all = fisher.test(m)$p.value, 
        pNaive_pBoth = fisher.test(m[, c(1, 2)])$p.value,
        mNaive_mBoth = fisher.test(m[, c(4, 5)])$p.value,
        pNaive_mMeta = fisher.test(m[, c(4, 6)])$p.value,
        pNavie_mNaive = fisher.test(m[, c(1, 4)])$p.value,
        pBoth_mBoth = fisher.test(m[, c(2, 5)])$p.value,
        pMeta_mMeta = fisher.test(m[, c(3, 6)])$p.value
        )
    })

fisher_result = data.frame(t(fisher_result))
dat_combind_sub = cbind(dat_combind_sub, fisher_result)

mut_Freq_test_result = dat_combind_sub
write_tsv(mut_Freq_test_result, "../data/03_mut_Freq_test_therapy.tsv")


#' # Draw graph
#+ mutation Freq by therapy,  fig.height=7, fig.width=15, dev='pdf'
dat_combind_long = dat_combind_sub[, .(symbol, freq_p_naive, freq_p_meta, freq_p_both, freq_m_naive, freq_m_meta, freq_m_both)] %>% melt
dat_combind_long %<>% rename(c("variable" = "group", "value" = "mutation_freq"))
dat_combind_long$symbol %<>% factor(levels = dat_combind_long[, .(s = sum(mutation_freq)), by=symbol][order(s, decreasing=T), symbol])

colors = brewer.pal(5, "Blues")[c(2,3,5)]
colors[4:6] =  brewer.pal(5, "Greens")[c(2,3,5)]
ggplot(data = dat_combind_long) +
    aes(x=symbol, y=mutation_freq, label=paste0(round(mutation_freq, 1), "%"), fill=group) +
    coord_cartesian(ylim = c(0, max(dat_combind_long$mutation_freq) * 1.1)) +
    geom_bar(data = dat_combind_long, position = "dodge", stat = "identity") +
    geom_text(position=position_dodge(width=.9), vjust=0.7, hjust=-0.2, angle=90, size=5) +
    scale_fill_manual(values=colors) + theme
