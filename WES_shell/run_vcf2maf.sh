#!/bin/bash

# config 
source "../config.sh"
conda activate wes-variants-calling
cd ${project_dir}

vcf2maf=/home/data/dataset/CHOL_2sample_new/bioinformatics-tools/vcf2maf-1.6.21.20230511/vcf2maf.pl
input_vcf="variants_calling/all_samples_2callers.vcf"
output_vcf="variants_calling/all_samples_2callers.vep${cache_version}1.vcf"
cache_version=111  
dir_cache=/home/data/database/vep

${vcf2maf} --input-vcf ${input_vcf} \
    --output-maf ${output_vcf} \
    --ref-fasta ${ref_fa} \
    --vep-data ${dir_cache} \
    --vep-path "/home/${USER}/miniconda3/envs/wes-variants-calling/share/ensembl-vep-111.0-0" \
    --cache-version ${cache_version} --vep-overwrite \
    --ncbi-build GRCh38 