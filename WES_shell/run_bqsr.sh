#!/bin/bash

bqsr(){
    # input: re.bam
    # intermiade: .re.sorted.bam, re.sorted.du.bam, re.sorted.du.matrix.txt, recalibration_report.txt
    # output: re.sorted.du.bqsr.bam
    sample_ID=${1}

    # step1. bam file排序
    java -Xmx50g -jar ${picard} SortSam \
		I="${input_bam_dir}/${sample_ID}.bam" \
		O="${output_bam_dir}/${sample_ID}.sorted.bam" \
		SO=coordinate \
        VALIDATION_STRINGENCY=STRICT \
        CREATE_INDEX=true

    # step2. 將重複序列進行標注
	java -jar ${picard} MarkDuplicates \
		I="${output_bam_dir}/${sample_ID}.sorted.bam" \
		O="${output_bam_dir}/${sample_ID}.sorted.du.bam" \
		M="${output_bam_dir}/${sample_ID}.sorted.du.matrix.txt" \
		REMOVE_DUPLICATES=false \
        CREATE_INDEX=true 2>> MARK_DUP.log

    # step3. 校正base品質
	gatk BaseRecalibrator \
		-R ${ref_fa} \
		-I "${output_bam_dir}/${sample_ID}.sorted.du.bam" \
		-O "${output_bam_dir}/${sample_ID}.recalibration_report.txt" \
		--known-sites ${dbsnp}

    # step4. 套用校正
    gatk ApplyBQSR \
		-R ${ref_fa} \
		-I "${output_bam_dir}/${sample_ID}.sorted.du.bam" \
		-O "${output_bam_dir}/${sample_ID}.sorted.du.bqsr.bam" \
		--bqsr-recal-file "${output_bam_dir}/${sample_ID}.recalibration_report.txt"
}

# config 
source "../config.sh"
conda activate wes-preprocessing
picard="/home/${USER}/miniconda3/envs/wes-preprocessing/share/picard-*/picard.jar" # 抓取使用者的環境資料夾

# project config
cd ${project_dir}
input_bam_dir="${project_dir}/bam"
output_bam_dir="${project_dir}/bqsr_bam"
mkdir -p ${output_bam_dir}
mkdir -p ${output_bam_dir}/log

case_IDs_file="${project_dir}/case_IDs.tsv"
while IFS= read -r case_ID
	do
		echo "Processing ${case_ID}"
		bqsr "${case_ID}-T" > "${output_bam_dir}/log/${case_ID}-T.bqsr.log" 2>&1 || { echo "Error processing ${case_ID}-T"; exit 1; }
		bqsr "${case_ID}-N" > "${output_bam_dir}/log/${case_ID}-N.bqsr.log" 2>&1 || { echo "Error processing ${case_ID}-N"; exit 1; }
		python "${tools_dir}/others/send_msg.py" "BQSR ${case_IDs} done"
	done < "$case_IDs_file"