#!/bin/bash
# config 
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate vep

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}

input_vcf="variants_calling/all_samples_2callers.vcf"
cache_version=103
output_vcf="variants_calling/all_samples_2callers.vep${cache_version}1.vcf"
dir_cache=/home/data/database/vep
assembly=GRCh38
threads=20
assembly_fasta=/home/data/data_thousand/gatk_index/Homo_sapiens_assembly38.fasta

vep -i $input_vcf \
    -o $output_vcf \
    --vcf \
    --cache --dir_cache $dir_cache --cache_version $cache_version \
    --assembly $assembly --force_overwrite \
    --fork $threads \
    --offline --no_stats --everything -pick \
    --fasta $assembly_fasta

