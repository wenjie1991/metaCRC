# input     ../data/02_Primary_mut_rate_right.tsv
#           ../data/02_Primary_mut_rate_left.tsv
#           ../data/02_Metastasis_mut_rate_right.tsv
#           ../data/02_Metastasis_mut_rate_left.tsv
# Output    FigureS2

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
dat_trios_p_right = fread("../data/02_Primary_mut_rate_right.tsv")
dat_trios_p_left = fread("../data/02_Primary_mut_rate_left.tsv")
dat_trios_m_right = fread("../data/02_Metastasis_mut_rate_right.tsv")
dat_trios_m_left = fread("../data/02_Metastasis_mut_rate_left.tsv")

n_right_trios = dat_trios_p_right[1, freq_n / freq * 100 %>% round]
n_left_trios = dat_trios_p_left[1, freq_n / freq * 100 %>% round]
n_right_trios
n_left_trios

#' # Combind data 
dat_combind = merge(dat_trios_p_right, dat_trios_p_left, all=T, by="symbol", suffixes=c("_p_right", "_p_left")) %>% 
    merge(dat_trios_m_right, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_m_right", "freq"="freq_m_right")) %>%
    merge(dat_trios_m_left, all=T, by="symbol") %>% rename(c("freq_n"="freq_n_m_left", "freq"="freq_m_left"))
dat_combind[dat_combind %>% is.na] = 0

#' # Output Table
dat_combind[freq_p_right > 5 | freq_p_left > 5 | freq_m_right > 5| freq_m_left > 5, .(symbol, freq_p_right, freq_p_left, freq_m_right, freq_m_left)]
write_tsv(dat_combind, "../data/03_mutation_frequency.tsv")

#' # Statistics test
dat_combind_sub = dat_combind[freq_p_right > 10 | freq_p_left > 10 | freq_m_right > 10| freq_m_left > 10]
fisher_result = sapply(1:nrow(dat_combind_sub), function(i) {
    m = dat_combind_sub[i, .(
        freq_n_p_right, n_right_trios - freq_n_p_right, 
        freq_n_p_left, n_left_trios - freq_n_p_left,
        freq_n_m_right, n_right_trios - freq_n_m_right,
        freq_n_m_left, n_left_trios - freq_n_m_left)] %>% unlist %>% matrix(nrow=2)
    c(
        all = fisher.test(m)$p.value, 
        pR_pL = fisher.test(m[, c(1, 2)])$p.value,
        mR_mL = fisher.test(m[, c(3, 4)])$p.value,
        pR_mR = fisher.test(m[, c(1, 3)])$p.value,
        pL_mL = fisher.test(m[, c(2, 4)])$p.value
        )
    })

fisher_result = data.frame(t(fisher_result))
dat_combind_sub = cbind(dat_combind_sub, fisher_result)

mut_Freq_test_result = dat_combind_sub
write_tsv(mut_Freq_test_result, "../data/03_mut_Freq_test_side.tsv")


#' # Draw graph
#+ mutation Freq by location,  fig.height=7, fig.width=15, dev='pdf'
dat_combind_long = dat_combind_sub[, .(symbol, freq_p_right, freq_p_left, freq_m_right, freq_m_left)] %>% melt
dat_combind_long %<>% rename(c("variable" = "group", "value" = "mutation_freq"))
dat_combind_long$symbol %<>% factor(levels = dat_combind_long[, .(s = sum(mutation_freq)), by=symbol][order(s, decreasing=T), symbol])

colors = brewer.pal(5, "Blues")[3:5]
colors[3:4] =  brewer.pal(4, "Greens")[c(2, 4)]
ggplot(data = dat_combind_long) +
    aes(x=symbol, y=mutation_freq, label=paste0(round(mutation_freq, 1), "%"), fill=group) +
    coord_cartesian(ylim = c(0, max(dat_combind_long$mutation_freq) * 1.1)) +
    geom_bar(data = dat_combind_long, position = "dodge", stat = "identity") +
    geom_text(position=position_dodge(width=.9), vjust=0.7, hjust=-0.2, angle=90, size=4) +
    scale_fill_manual(values=colors) + theme


