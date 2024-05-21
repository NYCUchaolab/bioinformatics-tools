#!/bin/bash

# function 


# config 
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate pindel

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
samplesheet_path="${project_dir}/gdc_sample_sheet.2024-02-26.tsv"

# package config and mkdir
pindel_dir="variants_calling/pindel2"

case_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))
# for sampleID in sample_IDs
for case_ID in "${case_IDs[@]}"
do
    echo "${case_ID}"
    tumor="${case_ID}_01A"
    normal="${case_ID}_11A"
    echo ${normal}
    echo ${tumor}
    # run_pindel ${case_ID}
done