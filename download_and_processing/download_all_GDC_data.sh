#!/bin/bash

# conda 
env_name="wes" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

# config 
project_dir="/home/data/dataset/LUAD_asian_8_sample" # 專案資料夾
source "${project_dir}/config.sh"
cd ${project_dir}

# mkdir
download_dir="${project_dir}/bam"
mkdir -p "${download_dir}"
mkdir -p "${download_dir}/tmp"
mkdir -p "${download_dir}/log"

# load sample sheet 
sample_sheet="${project_dir}/10_sample_sheet.tsv"
GDC_IDs=($(tail -n +2 ${sample_sheet} | cut -f 1 | sort -u))

# download all files to tmp
gdc-client download ${GDC_IDs[*]} \
    -n 20 \
    -d "${download_dir}/tmp" \
    -t ${GDC_token}
    --log-file "${download_dir}/log/download.log" \
    --latest 

# move and rename
tail -n +2 $sample_sheet | while IFS=$'\t' read -r -a sample
do
    GDC_ID=${sample[0]}
    file_name=${sample[1]}
    sample_ID=${sample[8]}
    mv ${download_dir}/tmp/${GDC_ID}/${file_name} ${download_dir}/${sample_ID}.bam
	mv ${download_dir}/tmp/${GDC_ID}/${file_name}.bai ${download_dir}/${sample_ID}.bam.bai
done 

# delete files
rm -r "${download_dir}/tmp