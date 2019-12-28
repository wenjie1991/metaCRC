# input         data/depth50/5000/
# output        Figure3D, Figure3E
#               ../data/VAP/result_depth50/result_dt.tsv               

library(stringr)
library(data.table)
# Load VAP
source("../VAP/playground/MRS_analysis_virtual_sample.R")

files = dir("data/depth50/5000/", full=T)


draw = function(in_file) {
    sampAB = fread(in_file) %>% as.data.frame
    snr = basename(in_file)
    samples = c(1:8)
    sampAB_new = pubOrSub.simu(sampAB, samples, minAF=0.05, minDepTotal=5*length(samples), groupName = "", pAF=0.25)
    depths = paste0("depth", samples)
    samples = paste0("maf", samples)
    plotRes.simVAF.matrix.pdf(sampAB_new, samples, depths, pdfsize = 40, plotType = "AF", snr=snr, sns=samples)
}

parameters = list()
for (in_file in files) {
    parameters[[in_file]] = draw(in_file)
}

par_names = names(parameters)
names(parameters[[1]])[1]

FST_v = c()
KSD_v = c()
result_dt = data.table()

for (l_i in names(parameters)) {
    for (l_j in names(parameters[[l_i]])) {
        matched = str_match(l_j, "simulMRS_deme(\\d+)_s(\\d+)percent_(\\d+)samples_u(\\d+\\.\\d+)_(\\d+)meta_init(\\d+)meta_cell_n(\\d+)depth(\\d+)\\.txt_(maf\\d+)_(maf\\d+)_stats")[1, -1]
        result_dt = rbind(result_dt, data.table(append(matched, parameters[[l_i]][[l_j]][c("FST", "KSD")] %>% unlist) %>% t))
    }
}

names(result_dt) = c("deme_n", "selection", "sample_n", "mutaton_rate", "batch", "meta_init", "meta_cell_n", "depth", "subA", "subB", "FST", "KSD")

write_tsv(result_dt, "../data/VAP/result_depth50/result_dt.tsv")
