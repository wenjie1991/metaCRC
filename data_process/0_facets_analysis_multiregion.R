# input         ../data/multiregion_biopsies/*
# output        ../data/multiregion_biopsies_facet_purity.tsv

patient_name_v = c("CHET54", "Ecol")

get_cval = function(patient_name) {
    c(CHET54 = 300, Ecol = 700)[patient_name]
}

estimate_purity = function(patient_name) {
    library(stringr) 
    ann_dt = fread(str_glue("../data/multiregion_biopsies/03_annotation/{patient_name}_local.hg19_multianno.csv"))
    af_dt = fread(str_glue("../data/multiregion_biopsies/03_annotation/{patient_name}.tsv"))

    normal_name = names(af_dt)[-(1:5)] %>% sub("\\.\\w+", "", .) %>% unique %>% grep("ND", ., invert=F, value=T)
    tumor_name_v = names(af_dt)[-(1:5)] %>% sub("\\.\\w+", "", .) %>% unique %>% grep("ND", ., invert=T, value=T)

    purity_v = c()

    do_fit = function(tumor_name) {
        normal_ref_sample = paste0(normal_name, ".REF")
        normal_alt_sample = paste0(normal_name, ".ALT")

        tumor_ref_sample = paste0(tumor_name, ".REF")
        tumor_alt_sample = paste0(tumor_name, ".ALT")


        af_dt_sub = af_dt[, c("CHROM", "POS_start", "POS_end", "REF", "ALT", normal_ref_sample, normal_alt_sample, tumor_ref_sample, tumor_alt_sample), with=F]
        names(af_dt_sub)[6:9] = c("normal_ref", "normal_alt", "tumor_ref", "tumor_alt")

        ann_dt_sub = ann_dt[Start == End, .(Chr, Start, End, Ref, Alt, genome1000 = `1000g2015aug_all`)]
        ann_dt_sub = ann_dt_sub[genome1000 > 0]

        setkeyv(af_dt_sub, c("CHROM", "POS_start", "POS_end"))
        af_dt_sub = af_dt_sub[ann_dt_sub[, 1:5, with=F]]

        af_dt_sub[, loc := paste0(CHROM, "_", POS_start, "_", POS_end)]
        dup_loc = af_dt_sub$loc[duplicated(af_dt_sub$loc)]
        af_dt_sub = af_dt_sub[!(loc %in% dup_loc)][normal_ref + normal_alt > 20 & tumor_ref + tumor_alt > 20]

        loc_dt = af_dt_sub[, .(Chromosome = CHROM %>% sub("chr", "", .), Position = POS_start)]

        df = data.frame(loc_dt, af_dt_sub[, .(
                NOR.DP = normal_alt + normal_ref,
                NOR.RD = normal_ref,
                TUM.DP = tumor_ref + tumor_alt,
                TUM.RD = tumor_ref
                )])
        library(facets)
        set.seed(1)
        print(tumor_name)
        xx = preProcSample(df, gbuild = "hg19")
        cval = get_cval(patient_name)
        oo = procSample(xx, cval = cval)
        fit = emcncf(oo, trace = T)
        list(x = oo, emfit = fit)
    }

    #     res = do_fit(tumor_name_v[3])
    #     plotSample(res$x, res$emfit)
    #     logRlogORspider(res$x$out, res$x$dipLogR)

    sapply(tumor_name_v, function(x) {
        do_fit(x)$emfit$purity
    }) -> purity_v
    names(purity_v) = tumor_name_v
    purity_v
    df = data.frame(
        sample = c(tumor_name_v, normal_name), 
        purity = c(purity_v, 0)
    )
    df
}

lapply(patient_name_v, estimate_purity) %>% rbindlist -> purity_df

write_tsv(purity_df, "../data/multiregion_biopsies_facet_purity.tsv")


