#!/bin/bash -

# input         ../data/pyclone/tsv_FACETS_cnv/
# output        ../data/pyclone/result/

set -o nounset                                  # Treat unset variables as an error
personID="$1"
tumor_purity_p="$2"
tumor_purity_m="$3"
mkdir -p result
mkdir -p result/
PyClone run_analysis_pipeline --in_files ../data/pyclone/tsv_FACETS_cnv/${personID}_primary.tsv ../data/pyclone/tsv_FACETS_cnv/${personID}_metastasis.tsv --working_dir ../data/pyclone/result/${personID} --prior total_copy_number --tumour_contents ${tumor_purity_p} ${tumor_purity_m}

