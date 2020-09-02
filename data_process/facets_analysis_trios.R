# input         ../data/matched_triols_annovar_output/*
# output        ../data/tumor_purity_facet.tsv

patient_name_v = dir("../data/matched_triols_annovar_output/02_annovar/") %>% sub("_.*", "", .) %>% unique

patient_name_v = c(patient_name_v %>% grep("AMC", ., value = T), patient_name_v %>% grep("ID\\d{7}", ., value = T))

get_cval = function() {
    150
}

# patient_name = patient_name_v[1]

do_fit = function(patient_name, tissue = "primary") {
    library(stringr) 
    patient_name2 = patient_name %>% sub("$-", "", .)
    ann_dt = fread(str_glue("../data/matched_triols_annovar_output/02_annovar/{patient_name}_local.hg19_multianno.csv"))
    af_dt = fread(str_glue("../data/matched_triols_annovar_output/03_genotype_field/{patient_name2}_genotype_field_AD_correct.tsv"))

    primary_ad = af_dt$C.AD %>% str_split_fixed(",", 2) %>% data.frame %>% data.matrix
    normal_ad = af_dt$N.AD %>% str_split_fixed(",", 2) %>% data.frame %>% data.matrix
    metastasis_ad = af_dt$M.AD %>% str_split_fixed(",", 2)  %>% data.frame %>% data.matrix

    if (tissue == "primary") {
        tumor_ad = primary_ad
    } else {
        tumor_ad = metastasis_ad
    }


    af_dt_sub      = af_dt[, .(CHROM, POS, POS, REF, ALT,
        norma_ref  = normal_ad[, 1],
        normal_alt = normal_ad[, 2],
        tumor_ref  = tumor_ad[, 1],
        tumor_alt  = tumor_ad[, 2] 
        )]
    names(af_dt_sub) = c("CHROM", "POS_start", "POS_end", "REF", "ALT", "normal_ref", "normal_alt", "tumor_ref", "tumor_alt")

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
    print(patient_name)
    xx = preProcSample(df, gbuild = "hg19")
    cval = get_cval()
    oo = procSample(xx, cval = cval)
    fit = emcncf(oo, trace = T)
    list(x = oo, emfit = fit)

}

library(stringr)
sapply(patient_name_v, function(x) {
    res = do_fit(x, "primary")
    patient_name_x = x %>% sub("-$", "", .)
    write_tsv(res$emfit$cncf, str_glue("../data/FACETS/cncf_{patient_name_x}_primary.tsv"))
    res$emfit$purity
}) -> purity_primary_v
sapply(patient_name_v, function(x) {
    res = do_fit(x, "metastasis")
    patient_name_x = x %>% sub("-$", "", .)
    write_tsv(res$emfit$cncf, str_glue("../data/FACETS/cncf_{patient_name_x}_metastasis.tsv"))
    res$emfit$purity
}) -> purity_metastasis_v

df = data.frame(
    personID = patient_name_v %>% sub("-$", "", .)
    , purity_p = purity_primary_v
    , purity_m = purity_metastasis_v
)
write_tsv(df, "../data/tumor_purity_facet.tsv")
