#!/bin/bash

# config 
source "../config.sh"
conda activate wes-variants-calling
cd ${project_dir}

input_vcf="variants_calling/all_samples_2callers.vcf"
output_vcf="variants_calling/all_samples_2callers.vep${cache_version}.vcf"
cache_version=103
dir_cache=/home/data/database/vep
assembly=GRCh38

vep -i ${input_vcf} \
    -o ${output_vcf} \
    --vcf \
    --cache --dir_cache ${dir_cache} --cache_version ${cache_version} \
    --assembly ${assembly} --force_overwrite \
    --fork ${threads} \
    --offline --no_stats --everything -pick \
    --fasta ${ref_fa}

