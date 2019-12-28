# input         ../data/VAP/
#               ../data/config.mapping.tsv
# output        Figure4A, Figure4B, Figure4C, Figure4D, Figure4E


library(gtools)
source("../VAP/analysis/MRS_analysis.R")

## load the samples
person_list = dir("../data/VAP/")
data_dir = "../data/VAP/"
cd = fread("../data/config.mapping.tsv")

## for each person
for (personID.i in person_list[25:37]) {
    # for (personID.i in person_list) {
    print(which(personID.i == person_list))
    print(personID.i)
    dir.create(personID.i)

    d1 = fread(paste0(data_dir, personID.i, "/mutect.snv.res.filtered.classified.founds.nopara.somatic.table"))
    #     d1 = d0[!((ref == "C" & alt == "T") | (ref == "G" & alt == "A"))]
    d.df = d1 %>% data.frame

    samples = cd[personID == personID.i, sampleID] %>% grep("sn", ., invert=T, value=T)
    if (personID.i == "CHET58") {
        samples = samples[-3]
    }
    normal = "sn"
    original = 4.3  # 4.3
    cmedianTh = 2  # 2

    res1 = getSampMutMulti(samples, normal, d.df, cmedianTh, original)
    res2 = adjust.ccf.titan.multi(res1, samples, 0.10, titanPath = paste0(data_dir, personID.i, "/titan/"))
    res3 = pubOrSub(res2, samples)

    cb = combinations(length(samples), 2, samples)
    for (i in 1:nrow(cb)) {
        plotRes.multi.pdf(res3, sampName=paste0(personID.i, "/", cb[i, 1], "_", cb[i, 2]), sn1 = paste0(cb[i, 1], "mafa"), sn2=paste0(cb[i, 2], "mafa"), minAF = 0.07)
    }
}
