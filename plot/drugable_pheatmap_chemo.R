# input         ../data/02_mutDT.tsv
#               ../data/phenotype.csv
# output        FigureS6

library(data.table)
library(magrittr)

tab = fread("../data/02_mutDT.tsv", key="sample")
cd = fread("../data/phenotype.csv", key="pathoID")
cd = cd[(metastasis_site %in% c("Liver"))]
cd[, personID := pathoID]
cd[grepl("^\\d+$", pathoID), personID := paste0("ID", substr(pathoID, 3, 9))]
drug_gene = fread("../data/allActionableVariants.txt")

tab = merge(tab, cd, all=T, by.x="sample", by.y="personID")

mut_heatmap = function(mutDT) {
    library(RColorBrewer)
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuBu")))(100)[100:1][c(1, 2:100)]

    match_color = function(x, range, color){
        color[ceiling((length(color) - 1) * (x - range[1]) / (range[2] - range[1])) + 1]
    }


    draw_mut = function(x, y, mut_width, mut_height, AF, color) {
        n = length(AF)
        loci_width = mut_width / n
        xleft = x + cumsum(c(rep(loci_width, n))) - loci_width
        xright = xleft + loci_width
        ybottom = rep(y, n)
        ytop = ybottom + mut_height
        for (i in 1:n) {
            col = match_color(AF[i], c(0, 1), color)
            rect(xleft[i], ybottom[i], xright[i], ytop[i], border = F, col = col)
        }
    }

    add_color_legend = function(x, y, breaks, width = 10) {
        for (i in 1:length(breaks)) {
            x1 = x
            y1 = y + i * 10 + 1
            draw_mut(x1, y1, mut_width = 10, mut_height = 10, AF = as.numeric(breaks[i]), color = color)
            text(x1 + width, y = y1 + width / 2, labels = breaks[i], pos = 4)
        }
    }

    draw_heatmap = function(mutDT) {
        # mutDT: sample, gene, loci, AF1, AF2

        samples = levels(mutDT$sample)
        genes = levels(mutDT$gene)

        mut_width = 12 
        mut_height = 12
        xRange = c(0, c(mut_width + 2) * length(samples) + 12 + 100)
        yRange = c(0, (mut_height + 2 + mut_height + 4) * length(genes) + 100)
        frame()
        plot.window(xRange, yRange)

        add_color_legend(x = xRange[2] - 20, y = yRange[2] /2, breaks = c("0", "0.25", "0.50", "0.75", "1.00"), width = 10)

        ## draw mutation
        for (s in 1:length(samples)) {
            for (g in 1:length(genes)) {
                x = (s - 1)* (mut_width + 2) + 20
                y1 = (g - 1) * (mut_height + 2 + mut_height + 4) + 40
                y2 = y1 - 2 - mut_height
                AF1 = mutDT[sample == samples[s] & gene == genes[g], AF1] %>% unlist
                AF2 = mutDT[sample == samples[s] & gene == genes[g], AF2] %>% unlist
                if (length(AF1) == 0) AF1 = 0
                if (length(AF2) == 0) AF2 = 0
                draw_mut(x, y1, mut_width, mut_height, AF = AF1, color)
                draw_mut(x, y2, mut_width, mut_height, AF = AF2, color)
            }
        }

        ## add patient annotation
        cd_sub = cd[, .(personID, Both_treated, primary_site)] %>% unique
        brewer.pal.info

        ### Therapy
        therapy = cd_sub$Both_treated
        therapy = factor(cd_sub$Both_treated, levels = c("Chemonaive", "Metastasis_treated", "Both_treated", NA), exclude=NULL)
        therapy_colors = c(brewer.pal(5, "Oranges")[c(1, 2, 4)], rgb(0, 0, 0, 0.7))[as.integer(therapy)]
        names(therapy_colors) = cd_sub$personID
        for (s in 1:length(samples)) {
            x = (s) * (mut_width + 2) + 20
            y1 = (length(genes) + 1) * (mut_height + 2 + mut_height + 4) + 40
            y2 = y1 - 2 - mut_height 
            rect(xleft=x-mut_width, ybottom=y1, xright=x, ytop=y2, col=therapy_colors[samples[s]])
        }

        ### Location 
        location = cd_sub$primary_site
        location = factor(location, levels=c("left", "right", NA), exclude=NULL)
        location_colors = c(brewer.pal(3, "Greens")[2], brewer.pal(3, "Blues")[2], rgb(0, 0, 0, 0.7))[as.integer(location)]
        names(location_colors) = cd_sub$personID
        for (s in 1:length(samples)) {
            x = (s) * (mut_width + 2) + 20
            y1 = (length(genes)) * (mut_height + 2 + mut_height + 4) + 40
            y2 = y1 - 2 - mut_height 
            rect(xleft=x-mut_width, ybottom=y1, xright=x, ytop=y2, col=location_colors[samples[s]])
        }
        

        ## add sample type
        colors = brewer.pal(5, "Greens")
        for (g in 1:length(genes)) {
            x1 = 6
            x2 = x1 + mut_width
            y1 = (g - 1) * (mut_height + 2 + mut_height + 4) + 40
            y2 = y1 - 2 - mut_height
            rect(xleft = x1, ybottom = y1, xright = x2, ytop = y1 + mut_height, col = colors[3])
            rect(xleft = x1, ybottom = y2, xright = x2, ytop = y2 + mut_height, col = colors[5])
        }

        
        ## add sample id
        for (s in 1:length(samples)) {
            x = s * (mut_width + 2) + 15
            y = 20
            text(x, y, labels = samples[s], srt = 60, pos = 2)
        }
        ## add gene symbol
        for (g in 1:length(genes)) {
            x = 0
            y1 = (g - 1) * (mut_height + 2 + mut_height + 4) + 40
            y2 = y1 - 2 - mut_height
            y = (y1 + y2) / 2 + mut_height / 2
            text(x, y, labels = genes[g], pos = 2)
        }

        ## add dash line to seperate gene
        for (g in 1:(length(genes) - 1)) {
            x1 = (0)* (mut_width + 2) + 20
            x2 = (length(samples))* (mut_width + 2) + 20
            y = (g - 1) * (mut_height + 2 + mut_height + 4) + 40 + 14
            lines(x = c(x1, x2), y = c(y, y), lty = 2)
        }
    }
    draw_heatmap(mutDT)
}


save_pdf= function(tab_sub, filename, gene_order, all_samples) {

    ### reorder
    m = tab_sub[, .(sample, gene)] %>% unique %>% dcast(gene ~ sample) %>% as.matrix
    rownames(m) = m[, 1]
    m = m[, -1]
    m[!is.na(m)] = 1
    m[is.na(m)] = 0 
    dimname = dimnames(m)
    m %<>% apply(2, as.numeric)
    dimnames(m) = dimname

    ### reorder gene 
    geneCount = apply(m, 1, sum)
    # geneCount["SYNE1"] = 30 + geneCount["SYNE1"]  ## change the SYNE1 rank
    m = m[order(geneCount), ]

    ## reorder mutation
    for (i in 1:nrow(m)) {
        m = m[, order(m[i, ])]
    }

    tab_sub$gene %<>% factor(levels = gene_order)
    sample_level = colnames(m) %>% rev
    sample_level = append(sample_level, all_samples[!(all_samples %in% sample_level)])
    tab_sub$sample %<>% factor(levels = sample_level)

    ## draw mutation heatmap
    pdf(filename, height = 20, width = 30)
    mut_heatmap(tab_sub)
    dev.off()
}

# genes = tab[, .(N = length(unique(sample))), by = gene][N >= 7, gene] %>% append("CTNNB1")
genes = tab[gene %in% drug_gene[Level %in% c("R1", "1", "2A", "2B"), Gene], gene] %>% unique

## Order matrix
tabNew = tab[gene %in% genes]
write_tsv(tabNew, "../data/05_pair_unconsist.tsv")
### reorder
m = tabNew[, .(sample, gene)] %>% unique %>% dcast(gene ~ sample) %>% as.matrix
rownames(m) = m[, 1]
m = m[, -1]
m[!is.na(m)] = 1
m[is.na(m)] = 0
dimname = dimnames(m)
m %<>% apply(2, as.numeric)
dimnames(m) = dimname

### reorder gene 
geneCount = apply(m, 1, sum)
# geneCount["SYNE1"] = 30 + geneCount["SYNE1"]  ## change the SYNE1 rank
m = m[order(geneCount), ]

## reorder mutation
for (i in 1:nrow(m)) {
    m = m[, order(m[i, ])]
}

gene_order = rownames(m) %>% na.omit


tabNew = tab[gene %in% genes & Both_treated=="Chemonaive"]
save_pdf(tabNew, "./drug_mut_heatmap_chemonaive.pdf", gene_order, cd[Both_treated == "Chemonaive", personID])

tabNew = tab[gene %in% genes & Both_treated=="Metastasis_treated"]
save_pdf(tabNew, "./drug_mut_heatmap_metastasis_treated.pdf", gene_order, cd[Both_treated == "Metastasis_treated", personID])

tabNew = tab[gene %in% genes & Both_treated=="Both_treated"]
save_pdf(tabNew, "./drug_mut_heatmap_both_treated.pdf", gene_order, cd[Both_treated == "Both_treated", personID])
