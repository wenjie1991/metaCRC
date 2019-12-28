#!/bin/bash -

# input         ../data/pyclone/tsv/
# output        ../data/pyclone/result/ 

set -o nounset                                  # Treat unset variables as an error
personID="$1"
mkdir -p result
mkdir -p result/
PyClone run_analysis_pipeline --in_files ../data/pyclone/tsv/${personID}_primary.tsv ../data/pyclone/tsv/${personID}_metastasis.tsv --working_dir ../data/pyclone/result/${personID} --prior total_copy_number

