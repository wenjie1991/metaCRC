prepare_Pyclone_input2.R
> input         ../data/01_prepare_data.allele_frequency.output.tsv
> output        ../data/pyclone/tsv/


run_pyclone.sh
> input         ../data/pyclone/tsv/

2_prior_mut_CC.R
> input         ../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv
>               ../data/CC_annovar_output/data_clinical_sample.txt
>               ../data/01_CC_somaticMutation_raw.txt
>               ../data/01_CC_somaticMutation.txt
> output
>               ../data/02_early_stage_CC_mut_rate.tsv
>               ../data/02_meta_stage_CC_mut_rate.tsv
>               ../data/02_meta_CC_mut_rate.tsv


1_CC_somatic_mutation.R
> Input         ../data/CC_annovar_output/03_genotype_field/
>               ../data/CC_annovar_output/02_annovar/
>               ../data/CC_annovar_output/02_annovar/
> Output        ../data/01_CC_somaticMutation_raw.txt
>               ../data/01_CC_somaticMutation.txt

prepare_Pyclone_input1.R
> input         ../data/01_trios_somaticMutation_raw.txt
> output        ../data/01_prepare_data.allele_frequency.output.tsv

2_prior_mut_trios.R
> input         ../data/01_trios_somaticMutation_raw.txt
>               ../data/01_trios_somaticMutation.txt
>               ../data/capture_kit_capture_region_files/MSK-IMPACT_panel.csv
> output
>               ../data/02_VAF_t_test.tsv
>               ../data/02_Primary_mut_rate.tsv
>               ../data/02_Metastasis_mut_rate.tsv
>               ../data/02_Primary_mut_rate_right.tsv
>               ../data/02_Primary_mut_rate_left.tsv
>               ../data/02_Metastasis_mut_rate_right.tsv
>               ../data/02_Metastasis_mut_rate_left.tsv
>               ../data/02_Primary_mut_rate_syn.tsv
>               ../data/02_Primary_mut_rate_sub.tsv
>               ../data/02_Metastasis_mut_rate_syn.tsv
>               ../data/02_Metastasis_mut_rate_sub.tsv
>               ../data/02_Primary_mut_rate_both.tsv
>               ../data/02_Primary_mut_rate_naive.tsv
>               ../data/02_Primary_mut_rate_meta.tsv
>               ../data/02_Metastasis_mut_rate_both.tsv
>               ../data/02_Metastasis_mut_rate_naive.tsv
>               ../data/02_Metastasis_mut_rate_meta.tsv
>               ../data/02_Primary_mut_rate_MSK.tsv
>               ../data/02_Primary_mut_rate_Local.tsv
>               ../data/02_Primary_mut_rate_Oncotarget.tsv
>               ../data/02_Metastasis_mut_rate_MSK.tsv
>               ../data/02_Metastasis_mut_rate_Local.tsv
>               ../data/02_Metastasis_mut_rate_Oncotarget.tsv
>               ../data/02_mutDT.tsv


1_matched_trios_somatic_mutation.R
> input     ../data/matched_triols_annovar_output/03_genotype_field/
>           ../data/matched_triols_annovar_output/02_annovar/
>           ../data/phenotype.csv
> output    ../data/01_trios_somaticMutation_raw.txt
>           ../data/01_trios_somaticMutation.txt
>           ../data/01_tumor_purity.tsv


0_analysis_coverage.R
> Input     ../data/VAP_QC/coverage*
> Output    ../data/00_coverage_greater_than_n.tsv
>           ../data/00_mean_coverage.tsv


