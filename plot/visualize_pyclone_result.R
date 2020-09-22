# input         ../data/pyclone/
# output        FigureS9

library(Cairo)
library(data.table)
library(magrittr)
library(ggplot2)
library(readr)

theme0 <- theme_bw() + theme(
    text = element_text(size = 15),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 0, vjust = 1, size = 15, hjust = 0.5),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    plot.title = element_text(hjust = 0.5),
    legend.position = "center",
) 

personList = "../data/pyclone//ID1116012 ../data/pyclone//ID1117428 ../data/pyclone//ID1126359 ../data/pyclone//AMC6 ../data/pyclone//AMC7 ../data/pyclone//AMC11 ../data/pyclone//AMC12 ../data/pyclone//AMC16 ../data/pyclone//AMC18 ../data/pyclone//AMC19 ../data/pyclone//AMC23 ../data/pyclone//ID1115829 ../data/pyclone//ID1119793 ../data/pyclone//ID1128834 ../data/pyclone//ID1135885" %>% strsplit(" ") %>% extract2(1) %>% basename



#+ output, fig.width=5, fig.height=4, dev='svg'
for (i in 1:length(personList)) {
	personID = personList[i]

	clone_cluster_table = fread(paste0("../data/pyclone/result/", personID, "/tables/cluster.tsv"))
	sample_information = clone_cluster_table$sample_id %>% strsplit("_") %>% ldply
	clone_cluster_table$personID = sample_information[[1]]
	clone_cluster_table$sampleType = sample_information[[2]]

	clone_cluster_table$cluster_id %<>% as.factor
	clone_cluster_table$sampleType %<>% revalue(., replace = c("primary" = "Primary", "metastasis" = "Metastasis"))

	p = ggplot(clone_cluster_table[size >= 5], aes(x = sampleType, color = cluster_id)) + 
		geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width=0.2, alpha=0.7) +
		geom_point(aes(y = mean), size=5, alpha=0.5) + 
		geom_line(aes(y = mean, group = cluster_id), alpha=0.7) + 
		labs(color = "Cluster ID", title = personID) +
		xlab("Sample Type") + ylab("CCF") + theme0
	print(p)
}
