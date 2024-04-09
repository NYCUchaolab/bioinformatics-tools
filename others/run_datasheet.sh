#!/bin/bash

# function 
# function (){}

# config 
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate wes

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
samplesheet_path="${project_dir}/gdc_sample_sheet.2024-02-26.tsv"

# package config and mkdir
# muse_dir="variants_calling/MuSE"
# mkdir -p variants_calling
# mkdir -p variants_calling/MuSE

# col 6 sample_ID TCGA-4G-AAZR sort 並去除重複
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# for sampleID in sample_IDs
for sample_ID in "${sample_IDs[@]}"
do
    echo "${sample_ID}"
    tumor="${sample_ID}_01A"
    normal="${sample_ID}_11A"
    echo ${normal}
    echo ${tumor}
done