# input         ../data/depth50_draw/
# output        Figure3B, Figure3C

library(stringr)
library(RColorBrewer)
library(Matrix)

colors = c("grey", RColorBrewer::brewer.pal(6, "YlGnBu"))
cluster_used_n = length(colors)

init = c("10000", "100000", "1000000", "10000000")
cell_n = c("1", "10", "100", "1000")
batch = c("0", "1")

draw_tumor = function(init_x, cell_n_x, batch_x)  {
    #     init_x = init[1]
    #     cell_n_x = cell_n[3]
    #     batch_x = batch[1]
    #     init_x = "1000000"
    #     cell_n_x = "1000"
    #     batch_x = "0"

    primary_dat = fread(paste0("../data/depth50_draw/CloneMap3D_peri_u5000_s0percent_8samples_u0.6_", batch_x, "meta_init", init_x, "meta_cell_n", cell_n_x, "depth50_primary.txt"))
    meta_dat = fread(paste0("../data/depth50_draw/CloneMap3D_peri_u5000_s0percent_8samples_u0.6_", batch_x, "meta_init", init_x, "meta_cell_n", cell_n_x, "depth50_metastasis_.txt"))

    primary_dat[, type := "primary"]
    meta_dat[, type := "meta"]
    dat = rbind(primary_dat, meta_dat)

    #     primary_dat[, last := str_split(lineage, "-") %>% (function(x) {sapply(x, function(y) {y[length(y)]})})]
    #     meta_dat[, last := str_split(lineage, "-") %>% (function(x) {sapply(x, function(y) {y[length(y)]})})]

    lineage = str_split(dat$lineage, "-")
    lineage_n = lineage %>% sapply(length)
    i = rep(1:length(lineage), lineage_n)
    j = as.integer(ordered(as.integer(unlist(lineage))))
    A = sparseMatrix(i, j, x=1)
    B = A[, colSums(A) > 1]
    set.seed(0)
    cluster_n = 8
    K = kmeans(B, cluster_n)

    clone_cluster = K$cluster
    cluster_tab = table(clone_cluster)
    #     lineage[clone_cluster == 2]
    cluster_introsize_median = tapply(lineage_n, clone_cluster, median)
    minor = names(cluster_tab)[order(cluster_tab)[1:(cluster_n - cluster_used_n)]]
    color_map = colors
    major_cluster_size_median = cluster_introsize_median[setdiff(names(cluster_introsize_median), minor)]
    names(color_map) = names(major_cluster_size_median)[order(major_cluster_size_median)]
    dat[, clone_cluster := as.character(clone_cluster)]
    dat[, clone_color := "white"]

    # intersect(unique(meta_dat$last), unique(primary_dat$last))
    dat[, clone_color := color_map[clone_cluster]]
    #     table(dat$clone_color)

    file_name = paste0("primary_metastasis_", init_x, "init_", cell_n_x, "meta_cell_n_", batch_x, ".svg")
    svg(paste0("./result_depth50_draw/", file_name), height=7/2*1.2, width=12/2)
    par(bty='n')
    with(dat[type == "primary"], plot(x = x - 60, y = y - 60, col = clone_color, xlim = c(-50, 160), ylim = c(-50, 50), pch = 15, cex=0.25, xaxt='n', yaxt='n', xlab="", ylab=""))
    with(dat[type == "meta"], points(x = x - 60 + 80, y = y - 60, col = clone_color, pch = 15, cex=0.25))
    axis(1, at = c(0, 80), labels = c("Primary", "Metastasis"), cex.axis=1.5)
    #     axis(2, cex.axis=1.5)
    xl = 130; ybottom = -20; ytop = 20;
    yb = seq(ybottom, ytop, length = cluster_used_n)
    yt = yb + (ytop - ybottom) / cluster_used_n * 0.8
    xr = xl + 5
    rect( xl, yb, xr, yt, col=colors)
    text(x = xr + 2, y = yb+2, labels = sort(major_cluster_size_median), adj = 0)
    text(x = xl, y = ytop+10, labels = "Mutation \nNunmber:\n", adj=0, cex=1.1)
    dev.off()
}

for (init_x in init) {
    for (cell_n_x in cell_n) {
        for (batch_x in batch) {
            draw_tumor(init_x, cell_n_x, batch_x)
        }
    }
}

