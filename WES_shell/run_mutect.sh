#!/bin/bash

# config 
source "../config.sh"
source "variants_calling"
conda activate wes-variants-calling
cd ${project_dir}

# project config
bam_dir="${project_dir}/bqsr_bam"

# package config and mkdir
mutect_dir="variants_calling/mutect"
mkdir -p variants_calling && chmod a+rw variants_calling
mkdir -p variants_calling/mutect && chmod a+rw variants_calling/mutect
mkdir -p variants_calling/mutect/log && chmod a+rw variants_calling/mutect/log

case_IDs_file="${project_dir}/case_IDs.tsv"
while IFS= read -r case_ID
	do
		echo "Processing ${case_ID}"
        run_mutect ${case_ID} > "${mutect_dir}/log/${case_ID}-T.mutect.log" 2>&1 || { echo "Error processing ${case_ID}"; exit 1; }
		python "${tools_dir}/others/send_msg.py" "mutect ${case_ID} done"
	done < "$case_IDs_file"