# input         ../data/01_trios_somaticMutation_raw.txt
# output        Figure S10A, Figure S10B, Figure S10C

tab = fread("../data/01_trios_somaticMutation_raw.txt")

tab$personID %>% unique 

personID_list = tab$personID %>% unique
personID_list %<>% grep("^ID", ., value=T)

# personID = personID_list[1]

theme0 <- theme_bw() + theme(
    text = element_text(size = 15),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 60, vjust = 1, size = 15, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "right",
) 

theme1 <- theme_bw() + theme(
    text = element_text(size = 15),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 0, vjust = 1, size = 15, hjust = 0.5),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "right",
) 

get_mut_type = function(personID_i, region_specific = F) {
    tab_sub = tab[personID == personID_i]
    print(tab_sub %>% nrow)

    if (region_specific == F) {
        is_primary_mut = tab_sub[, C.A / (C.A + C.R) > 0.25 & C.A + C.R > 20]
        is_metastasis_mut = tab_sub[, M.A / (M.A + M.R) > 0.25 & M.A + M.R > 20]
    } else {
        is_primary_mut = tab_sub[, C.A / (C.A + C.R) > 0.25 & C.A + C.R > 20 & M.A / (M.A + M.R) < 0.05 ]
        is_metastasis_mut = tab_sub[, M.A / (M.A + M.R) > 0.25 & M.A + M.R > 20 & C.A / (C.A + C.R) < 0.05]
    }

    mutation_mapping = c(
        "A/T" = "T/A",
        "A/G" = "T/C",
        "A/C" = "T/G",
        "C/A" = "G/T",
        "C/T" = "G/A",
        "C/G" = "G/C"
    )

    clean_mutation = function(is_mut) {
        x_mut_type = tab_sub[is_mut, paste0(REF, "/", ALT)]
        x_mut_type %<>% revalue(mutation_mapping)
        x_mut_type = x_mut_type[x_mut_type %in% mutation_mapping]
        x_mut_type
    }

    primary_mut_type_v = clean_mutation(is_primary_mut)
    metastasis_mut_type_v = clean_mutation(is_metastasis_mut)

    data.table(
        mut_type = c(primary_mut_type_v, metastasis_mut_type_v),
        sample_type = rep(c("Primary", "Metastasis"), c(length(primary_mut_type_v), length(metastasis_mut_type_v))),
        personID = personID_i
    )

}

## The mutation in primary and metastasis respectively
lapply(personID_list, function(x) {
    get_mut_type(x)
    }) %>% rbindlist -> mut_type_df


g = ggplot(mut_type_df) + aes(x = personID, fill = mut_type) 
p = g + geom_bar(position = "fill", stat = "count") + facet_grid(sample_type ~ .) + theme0
p + labs(x = "Patient ID", y = "Mutation type %", fill = "Mutation Type")

g = ggplot(mut_type_df) + aes(x = personID, fill = mut_type) 
p = g + geom_bar(position = "stack", stat = "count") + facet_grid(sample_type ~ .) + theme0
p + labs(x = "Patient ID", y = "Mutation type count", fill = "Mutation Type")



## The primary or metastasis specific mutation
lapply(personID_list, function(x) {
    get_mut_type(x, region_specific = T)
}) %>% rbindlist -> mut_type_df

d = mut_type_df
d$time = 16 - substr(d$personID, 3, 4) %>% as.numeric

d = d[personID != "ID0910513", .N, by = .(mut_type, time, personID)]
d$time %<>% as.factor


g = ggplot(mut_type_df) + aes(x = personID, fill = mut_type) 
p = g + geom_bar(position = "fill", stat = "count") + facet_grid(sample_type ~ .) + theme0
p + labs(x = "Patient ID", y = "Mutation type %", fill = "Mutation Type")

g = ggplot(d) + aes(x = time, y = N, fill = mut_type) 
# p = g + geom_bar(stat = "identity", position = "dodge", aes(group = mut_type)) + theme0
p = g + geom_boxplot(aes(x = time, y = N)) + theme1
p = p + labs(x = 'Sample "age" (year)', y = "Regional specific hVAF mutation Num", fill = "Mutation Type")
p

