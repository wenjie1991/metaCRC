# Input     ../data/VAP_QC/coverage*
# Output    ../data/00_coverage_greater_than_n.tsv
#           ../data/00_mean_coverage.tsv

library(data.table)
library(magrittr)
library(stringr)
library(plyr)

library(RColorBrewer)
colors = brewer.pal(5, "Blues")[c(3, 5)]
colors[3:5] =  brewer.pal(5, "Greens")[c(3:5)]

theme <- theme_bw() + theme(
    text = element_text(size = 25),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 0, vjust = 1, size = 15, hjust = 0.5),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "top",
    ) 

files <- dir(path = "../data/VAP_QC/coverage", pattern = "coverage", full.names = TRUE)  %>% grep("111999|0825371", ., invert=T, value=T)
labs <- str_match(basename(files), "(.*?)\\.")[, 2]

cummCoverage <- function(coverageFileData){
    (1 - cumsum(coverageFileData$V5))
}

covCumul <- lapply(files, function(f) {
    fread(f) %>% cummCoverage
})

## Mean covrage
summary_coverage = sapply(files, function(f) {
    fread(f) %>% (function(tab) {sum(tab$V2 * tab$V5)})()
}) 
summary_coverage = data.table(sample_name = labs, summary_coverage)
names(summary_coverage)[-1] = c("mean")
write_tsv(summary_coverage, "../data/00_mean_coverage.tsv")

summary_coverage[mean > 50]


# Coverage greater than 80
coverage_gt10 = sapply(files, function(f) {
    fread(f) %>% (function(tab) {1 - tab[V2 <= 10, sum(V5)]})()
}) 
coverage_gt20 = sapply(files, function(f) {
    fread(f) %>% (function(tab) {1 - tab[V2 <= 20, sum(V5)]})()
}) 
coverage_gt50 = sapply(files, function(f) {
    fread(f) %>% (function(tab) {1 - tab[V2 <= 50, sum(V5)]})()
}) 
coverage_gtn = data.table(sample_name = labs, coverage_gt10, coverage_gt20, coverage_gt50)
names(coverage_gtn)[2:4] = c("Coverage > 10", "Coverage > 20", "Coverage > 50")
coverage_long = melt(coverage_gtn)
write_tsv(coverage_long, "../data/00_coverage_greater_than_n.tsv")

