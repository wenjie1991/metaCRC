multiregion_VAF.R
```
# input         ../data/multiregion_biopsies/03_annotation/
#               ../data/multiregion_biopsies_facet_purity.tsv
#               ../data/multiregion_purbyes_estimateion.tsv
# output        Figure2C, Figure2D
```
pheatmap_chemo.R
```
# input         ../data/02_mutDT.tsv
#               ../data/phenotype.csv
# output        Figure1B
```
VAP_real_tumor_data_process.R
```
# input         ../data/VAP/
#               ../data/config.mapping.tsv
# output        Figure4A, Figure4B, Figure4C, Figure4D, Figure4E
```
vaf_plot.R
```
# input         ../data/01_trios_somaticMutation_raw.txt
#               ../data/01_trios_somaticMutation.txt")
# output        Figure2.A, Figure2.B
#               ../data/02_mutDT.tsv
```
draw_virtual_tumor.R
```
# input         ../data/depth50_draw/
# output        Figure3B, Figure3C
```
evaluate_pureBayes.R
```
## Legacy purity
## Purebayes with synonymous mutation as input
## Purebayes use FACETS CNV information
## Purebayes result use FACETS CNV infor & with synonymous mutation as input
## FACETS CNV
```
visualize_pyclone_result.R
```
# input         ../data/pyclone/
# output        FigureS9
```
FFPE_artifical_mutation.R
```
# input         ../data/01_trios_somaticMutation_raw.txt
# output        Figure S10A, Figure S10B, Figure S10C
```
mutation_frequency_therapy.R
```
# input     ../data/02_Primary_mut_rate_both.tsv
#           ../data/02_Primary_mut_rate_naive.tsv
#           ../data/02_Primary_mut_rate_meta.tsv
#           ../data/02_Metastasis_mut_rate_both.tsv
#           ../data/02_Metastasis_mut_rate_naive.tsv
#           ../data/02_Metastasis_mut_rate_meta.tsv
# output    FigureS4
```
mutation_frequency_stage.R
```
#input      ../data/02_early_stage_CC_mut_rate.tsv
#           ../data/02_meta_stage_CC_mut_rate.tsv
#           ../data/02_meta_CC_mut_rate.tsv
#           ../data/02_Primary_mut_rate.tsv
#           ../data/02_Metastasis_mut_rate.tsv
# output    Figure1A 
```
drugable_pheatmap_chemo.R
```
# input         ../data/02_mutDT.tsv
#               ../data/phenotype.csv
# output        FigureS6
```
VAP_real_tumor_statistics.R
```
# input             ../data/phenotype.csv, 
#                   ../data/00_mean_coverage.tsv
#                   ..data/patient_metrics.csv       
# output            Figure5A, Figure5B, Figure5C
```
VAP_virtual_tumor_data_process.R
```
# input         data/depth50/5000/
# output        Figure3D, Figure3E
#               ../data/VAP/result_depth50/result_dt.tsv               
```
tumor_purity_statistics.R
```
# input   ../data/01_tumor_purity.tsv
# output  Figure.S8 
```
mutation_frequency_location.R
```
# input     ../data/02_Primary_mut_rate_right.tsv
#           ../data/02_Primary_mut_rate_left.tsv
#           ../data/02_Metastasis_mut_rate_right.tsv
#           ../data/02_Metastasis_mut_rate_left.tsv
# Output    FigureS2
```
VAP_summary_statistics.R
```
# input         ../data/result_depth50/result_dt.tsv
#               ../data/result_depth50/result_dt.tsv
#               ../data/phenotype.csv
#               ../data/00_mean_coverage.tsv
#               ../data/VAP_real_data/patient_metrics.csv
# output        Figure5A, Figure5B, Figure5C
```
lilloplot.R
```
# input         ../data/01_trios_somaticMutation.txt
# output        FigureS3, FigureS7
```
mutation_frequency_resection_time.R
```
# input     ../data/02_Primary_mut_rate_syn.tsv
#           ../data/02_Primary_mut_rate_sub.tsv
#           ../data/02_Metastasis_mut_rate_syn.tsv
#           ../data/02_Metastasis_mut_rate_sub.tsv
# output    FigureS5
```
