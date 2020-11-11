run_pyclone.sh
```
#!/bin/bash -
```
ascat_analysis.R
```
# dt = fread("../data/multiregion_biopsies/03_annotation/CHET54_gene.hg19_multianno.csv")
```
matched_trios_purbayes.R
```
# Input         ../data/01_trios_somaticMutation_raw.txt
# Output        ../data/01_tumor_purity.tsv
#               ../data/01_tumor_purity_syn.tsv
#               ../data/01_tumor_purity_facetscnv.tsv
#               ../data/01_tumor_purity_facetscnvsyn.tsv
```
2_prior_mut_CC.R
```
# input         ../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv
#               ../data/CC_annovar_output/data_clinical_sample.txt
#               ../data/01_CC_somaticMutation_raw.txt
#               ../data/01_CC_somaticMutation.txt
# output
#               ../data/02_early_stage_CC_mut_rate.tsv
#               ../data/02_meta_stage_CC_mut_rate.tsv
#               ../data/02_meta_CC_mut_rate.tsv
```
multiregion_tumor_purity_purbyes.R
```
# input         ../data/multiregion_biopsies/03_annotation/
# output        ../data/multiregion_purbyes_estimateion.tsv
```
1_CC_somatic_mutation.R
```
# Input         ../data/CC_annovar_output/03_genotype_field/
#               ../data/CC_annovar_output/02_annovar/
#               ../data/CC_annovar_output/02_annovar/
# Output        ../data/01_CC_somaticMutation_raw.txt
#               ../data/01_CC_somaticMutation.txt
```
prepare_Pyclone_input1.R
```
# input         ../data/01_trios_somaticMutation_raw.txt
# output        ../data/01_prepare_data.allele_frequency.output2.tsv
```
2_prior_mut_trios.R
```
# input         ../data/01_trios_somaticMutation_raw.txt
#               ../data/01_trios_somaticMutation.txt
#               ../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv
# output
#               ../data/02_VAF_t_test.tsv
#               ../data/02_Primary_mut_rate.tsv
#               ../data/02_Metastasis_mut_rate.tsv
#               ../data/02_Primary_mut_rate_right.tsv
#               ../data/02_Primary_mut_rate_left.tsv
#               ../data/02_Metastasis_mut_rate_right.tsv
#               ../data/02_Metastasis_mut_rate_left.tsv
#               ../data/02_Primary_mut_rate_syn.tsv
#               ../data/02_Primary_mut_rate_sub.tsv
#               ../data/02_Metastasis_mut_rate_syn.tsv
#               ../data/02_Metastasis_mut_rate_sub.tsv
#               ../data/02_Primary_mut_rate_both.tsv
#               ../data/02_Primary_mut_rate_naive.tsv
#               ../data/02_Primary_mut_rate_meta.tsv
#               ../data/02_Metastasis_mut_rate_both.tsv
#               ../data/02_Metastasis_mut_rate_naive.tsv
#               ../data/02_Metastasis_mut_rate_meta.tsv
#               ../data/02_Primary_mut_rate_MSK.tsv
#               ../data/02_Primary_mut_rate_Local.tsv
#               ../data/02_Primary_mut_rate_Oncotarget.tsv
#               ../data/02_Metastasis_mut_rate_MSK.tsv
#               ../data/02_Metastasis_mut_rate_Local.tsv
#               ../data/02_Metastasis_mut_rate_Oncotarget.tsv
#               ../data/02_mutDT.tsv
```
facets_analysis_trios.R
```
# input         ../data/matched_triols_annovar_output/*
# output        ../data/tumor_purity_facet.tsv
```
supp_tumor_purity_qc.R
```
# tumor purity estimated by Purebayes,
```
1_matched_trios_somatic_mutation.R
```
# input     ../data/matched_triols_annovar_output/03_genotype_field/
#           ../data/matched_triols_annovar_output/02_annovar/
#           ../data/phenotype.csv
# output    ../data/01_trios_somaticMutation_raw.txt
#           ../data/01_trios_somaticMutation.txt
```
prepare_Pyclone_with_FACETS_cnv.R
```
# Input: ../data/01_prepare_data.allele_frequency.output2.tsv
# Output: ../data/pyclone/tsv_FACETS_cnv
```
facets_analysis_multiregion.R
```
# input         ../data/multiregion_biopsies/*
# output        ../data/multiregion_biopsies_facet_purity.tsv
```
facets_analysis.R
```
#     res = do_fit(tumor_name_v[3])
```
0_analysis_coverage.R
```
# Input     ../data/VAP_QC/coverage*
# Output    ../data/00_coverage_greater_than_n.tsv
#           ../data/00_mean_coverage.tsv
```
