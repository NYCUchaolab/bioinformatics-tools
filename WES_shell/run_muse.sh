#!/bin/bash

# config 
source "../config.sh"
source "variants_calling"
conda activate wes-variants-calling
cd ${project_dir}
bam_dir="${project_dir}/bqsr_bam"

# package config and mkdir
muse_dir="variants_calling/MuSE"
mkdir -p variants_calling && chmod a+rw variants_calling
mkdir -p variants_calling/MuSE && chmod a+rw variants_calling/MuSE
mkdir -p variants_calling/MuSE/log && chmod a+rw variants_calling/MuSE/log

case_IDs_file="${project_dir}/case_IDs.tsv"
while IFS= read -r case_ID
	do
		echo "Processing ${case_ID}"
        normal="${case_ID}-T"
        tumor="${case_ID}-N"
        cp -n "${bam_dir}/${normal}.sorted.du.bqsr.bai" "${bam_dir}/${normal}.sorted.du.bqsr.bam.bai"
        cp -n "${bam_dir}/${tumor}.sorted.du.bqsr.bai" "${bam_dir}/${tumor}.sorted.du.bqsr.bam.bai"

        run_muse "${case_ID}" > "${muse_dir}/log/${case_ID}-T.muse.log" 2>&1 || { echo "Error processing ${case_ID}"; exit 1; }
        python "${tools_dir}/others/send_msg.py" "mutect ${case_ID} done"
	done < "$case_IDs_file"