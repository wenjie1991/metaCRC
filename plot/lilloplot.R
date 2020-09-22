# input         ../data/01_trios_somaticMutation.txt
# output        FigureS3, FigureS7
library(trackViewer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(RColorBrewer)

lolliplot.plot = function(symbol, start_adj=0, end_adj=0, type='circle') {

    symbol = "APC"; end_adj = 10; start_adj = 10; type = 'circle'
    gene = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    id = select(org.Hs.eg.db, symbol, "ENTREZID", "SYMBOL")$ENTREZID
    gr = gene[id]

    mut = fread("../data/01_trios_somaticMutation.txt")
    mut_sub = mut[Gene.refGene == symbol]
    mut_sub[ExonicFunc.refGene == ".", ExonicFunc.refGene := Func.refGene]
    mut_sub[grepl("SNV", ExonicFunc.refGene), ExonicFunc.refGene := "nonsynonymous SNV"]
    mut_sub[grepl("frameshift|stopgain|splicing", ExonicFunc.refGene), ExonicFunc.refGene := "frameshift_indel/stopgain/splicing"]
    mut_sub$ExonicFunc.refGene %>% table 
    mut_sub_primary = mut_sub[C.A / (C.A + C.R + 1) > 0,
        .(Gene.refGene, CHROM, POS, REF, ALT, ExonicFunc.refGene, Both_treated, primary_site, Resection_timing, personID, value1 = C.A / (C.A + C.R + 1) * 100)]
    mut_sub_primary[, value2 := 100-value1]
    mut_sub_metastasis = mut_sub[M.A / (M.A + M.R + 1) > 0,  
        .(Gene.refGene, CHROM, POS, REF, ALT, ExonicFunc.refGene, Both_treated, primary_site, Resection_timing, personID, value1=M.A / (M.A + M.R + 1) * 100)]
    mut_sub_metastasis[, value2 := 100-value1]
    
    wilcox.test(mut_sub_primary[primary_site == "left", POS], mut_sub_primary[primary_site == "right", POS])
    mut_sub_primary[primary_site %in% c("left", "right")]

    mut_freq_primary = GRanges(mut_sub_primary$CHROM, IRanges(mut_sub_primary$POS, mut_sub_primary$POS))
    mut_freq_metastasis = GRanges(mut_sub_metastasis$CHROM, IRanges(mut_sub_metastasis$POS, mut_sub_metastasis$POS))

    mcols(mut_freq_primary) = DataFrame(
        REF                 = mut_sub_primary$REF
        , ALT               = mut_sub_primary$ALT
        , ExonicFunc.refGene= mut_sub_primary$ExonicFunc.refGene
        , Both_treated      = mut_sub_primary$Both_treated
        , primary_site      = mut_sub_primary$primary_site
        , Resection_timing  = mut_sub_primary$Resection_timing
        , score             = mut_sub_primary$value1
        #         , score             = 1
        , value1            = mut_sub_primary$value2
        , value2            = mut_sub_primary$value1
        #         , score             = mut_sub_primary$N
        )

    mcols(mut_freq_metastasis) = DataFrame(
        REF                    = mut_sub_metastasis$REF
        , ALT                  = mut_sub_metastasis$ALT
        , ExonicFunc.refGene   = mut_sub_metastasis$ExonicFunc.refGene
        , Both_treated         = mut_sub_metastasis$Both_treated
        , primary_site         = mut_sub_metastasis$primary_site
        , Resection_timing     = mut_sub_metastasis$Resection_timing
        , score                = mut_sub_metastasis$value1
        #         , score                 = 1
        , value1               = mut_sub_metastasis$value2
        , value2               = mut_sub_metastasis$value1
        #         , score                = mut_sub_metastasis$N
        )

    mut_freq_primary$SNPsideID="top"
    mut_freq_metastasis$SNPsideID="bottom"

    color_v = alpha(c("#000000", "#69b6ee"), c(0.7, 1))
    mut_freq_primary$color = color_v[factor(mut_freq_primary$ExonicFunc.refGene, levels=c("frameshift_indel/stopgain/splicing", "nonsynonymous SNV"))]
    mut_freq_metastasis$color = color_v[factor(mut_freq_metastasis$ExonicFunc.refGene, levels=c("frameshift_indel/stopgain/splicing", "nonsynonymous SNV"))]
    #     names(mut_freq_primar5) = mut_sub_primary$personID
    #     names(mut_freq_metastasis) = mut_sub_metastasis$personID

    mut_freq_both = c(mut_freq_primary, mut_freq_metastasis)


    mut_freq_both$border = alpha("#000000", 0.5)
    #     mut_freq_both$color = list(brewer.pal(5, "Blues")[c(3, 1)], brewer.pal(5, "Greens")[c(3,1)], brewer.pal(5, "Blues")[c(3, 1)], brewer.pal(5, "Blues")[c(3, 1)])[factor(mut_freq_both$Both_treated, levels = c("Chemonaive", "Metastasis_treated", "Both_treated", NA), exclude=NULL)]

    seqlevelsStyle(gr) = seqlevelsStyle(mut_freq_both) = "UCSC"

    trs = geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, gr=gr)
    print(names(trs)[1])

    x = trs[[1]]$dat
    x = x[x$feature == "CDS"]

    s = start(x)
    e = end(x)
    q_s = start(mut_freq_both)

    for (i in 1:length(q_s)) {
        to_start = q_s[i] - s
        to_end = e - q_s[i]
        is_in_exon = to_start > 0 & to_end > 0
        if (sum(is_in_exon) == 0) {
            q_s[i] = ifelse(min(abs(to_start)) > min(abs(to_end))
                , e[which.min(abs(to_end))]
                , s[which.min(abs(to_start))]
                )
        }
    }

    mut_freq_both@ranges = IRanges(start = q_s, end = q_s)

    mut_freq_tmp = pmapToTranscripts(mut_freq_both, GRangesList("tx_c" = x), ignore.strand=T)
    mcols(mut_freq_tmp) = mcols(mut_freq_both)
    mut_freq_both = mut_freq_tmp

    x = pmapToTranscripts(x, GRangesList("tx_c" = x), ignore.strand=T)

    if (as.character(strand(gr)[1]) == "-") {
        l = max(c(end(x), start(x)))
        mut_freq_both@ranges = IRanges(start=l-end(mut_freq_both), end=l-start(mut_freq_both))
        x@ranges = IRanges(start=l-end(x), end=l-start(x))
    }

    mut_freq_both@ranges = IRanges(start=start(mut_freq_both) / 3, end=end(mut_freq_both) / 3)
    x@ranges = IRanges(start=start(x)/3, end=end(x)/3)

    features = x
    names(features) <- c(trs[[1]]$name)
    features$fill <- "#d96918"
    #     features$height <- c(.03)

    
    #     r = range(mut_freq_both)
    r = range(x)
    end(r) = end(r) + end_adj
    start(r) = start(r) - start_adj
    #     r = gr
    #     r = range(mut_freq_both)

    ## Total needle plot
    print("########## Draw the needle plot for total samples ##########")
    mut_freq_location = mut_freq_both
    lolliplot(list(A=mut_freq_location), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5, dashline.col=NA)


    ## Location
    print("########## Primary Location ##########")
    mut_freq_location1 = mut_freq_both[mcols(mut_freq_both)$primary_site %in% "left"]
    mut_freq_location2 = mut_freq_both[mcols(mut_freq_both)$primary_site %in% "right"]
    lolliplot(list(A=mut_freq_location1), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5, dashline.col=NA)
    lolliplot(list(A=mut_freq_location2), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5, dashline.col=NA)
    #     dandelion.plot(mut_freq_location1, features, ranges=r, type=type, xaxis=T, cex=0.4)
    #     dandelion.plot(mut_freq_location2, features, ranges=r, type=type, xaxis=T, cex=0.4)
    # dandelion.plot(mut_freq, features, ranges=gr, type='pin', maxgaps=1/100)
        # statistics comparing the locaion between left and right

    ## Therapy
    #     print("########## Therapy ##########")
    mut_freq_therapy1 = mut_freq_both[mcols(mut_freq_both)$Both_treated %in% "Chemonaive"]
    mut_freq_therapy2 = mut_freq_both[mcols(mut_freq_both)$Both_treated %in% "Metastasis_treated"]
    mut_freq_therapy3 = mut_freq_both[mcols(mut_freq_both)$Both_treated %in% "Both_treated"]
    lolliplot(list(A=mut_freq_therapy1), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5)
    lolliplot(list(A=mut_freq_therapy2), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5)
    lolliplot(list(A=mut_freq_therapy3), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5)

    ## Resection timing
    #     print("########## Resection timing ##########")
    #     mut_freq_therapy1 = mut_freq_both[mcols(mut_freq_both)$Resection_timing %in% "Concurrent"]
    #     mut_freq_therapy2 = mut_freq_both[mcols(mut_freq_both)$Resection_timing %in% "Subsequent"]
    #     lolliplot(list(A=mut_freq_therapy1), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5)
    #     lolliplot(list(A=mut_freq_therapy2), list(x=features, y=features, z=features), ranges=r, type=type, xaxis=T, cex=0.5)
}


#' APC
#+ APC, fig.width=7, dev="pdf", fig.height=3
lolliplot.plot("APC", end_adj=10, start_adj=10)

#' TP53
# lolliplot.plot("TP53")
#+ TP53, fig.width=7, dev="pdf", fig.height=3
lolliplot.plot("TP53", end_adj=10, start_adj=10)

