#!/bin/bash

# config 
source "../config.sh"
source "variants_calling.sh"
conda activate wes-variants-calling
cd ${project_dir}
bam_dir="${project_dir}/bqsr_bam"

# package config and mkdir
pindel_dir="${project_dir}/variants_calling/pindel"
mkdir -p ${pindel_dir} && chmod a+rw ${pindel_dir}
mkdir -p ${pindel_dir}/log && chmod a+rw ${pindel_dir}/log
config_file="${pindel_dir}/${case_ID}.config.txt"

case_IDs_file="${project_dir}/case_IDs.tsv"
while IFS= read -r case_ID
do
echo "Processing ${case_ID}"
run_pindel ${case_ID} > "${pindel_dir}/log/${case_ID}.pindel.log" 2>&1 || { echo "Error processing ${case_ID}"; exit 1; }
python "${tools_dir}/others/send_msg.py" "pindel ${case_ID} done"
done < "$case_IDs_file"
