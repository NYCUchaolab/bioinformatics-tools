#!/bin/bash

run_bwa(){
    sample_ID=${1}
    
    bwa mem \
    -t ${threads} \
    -T 0 \
    -R "@RG\tID:${sample_ID}\tSM:${sample_ID}\tPL:ILLUMINA" \
    ${ref_fa} \
    "${fastq_dir}/${sample_ID}_R1.fq" \
    "${fastq_dir}/${sample_ID}_R2.fq" \
    | samtools view -Shb -o "${bam_dir}/${sample_ID}.bam" - 
}

# config 
source "../config.sh"
conda activate wes-preprocessing

cd ${project_dir}
fastq_dir="${project_dir}/fastq"
bam_dir="${project_dir}/bam"

case_IDs_file="${project_dir}/case_IDs.tsv"
while IFS= read -r case_ID
do
    echo "Processing ${case_ID}"
    run_bwa "${case_ID}-T" > "${bam_dir}/log/${case_ID}-T.bwa.log" 2>&1 || { echo "Error processing ${case_ID}-T"; exit 1; }
    run_bwa "${case_ID}-N" > "${bam_dir}/log/${case_ID}-N.bwa.log" 2>&1 || { echo "Error processing ${case_ID}-N"; exit 1; }
    python "${tools_dir}/others/send_msg.py" "bwa ${case_IDs} done"
done < "$case_IDs_file"