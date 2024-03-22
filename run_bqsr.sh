#!/bin/bash


bqsr(){
    # input: re.bam
    # intermiade: .re.sorted.bam, re.sorted.du.bam, re.sorted.du.matrix.txt, recalibration_report.txt
    # output: re.sorted.du.bqsr.bam

    sample_ID=${1}

    # step1. bam file排序
    java -Xmx50g -jar ${picard} SortSam \
		I="${input_bam_dir}/${sample_ID}.bam" \
		O="${output_bam_dir}/${sample_ID}.re.sorted.bam" \
		SO=coordinate \
        VALIDATION_STRINGENCY=STRICT \
        CREATE_INDEX=true

    # step2. 將重複序列進行標注
	java -jar ${picard} MarkDuplicates \
		I="${output_bam_dir}/${sample_ID}.re.sorted.bam" \
		O="${output_bam_dir}/${sample_ID}.re.sorted.du.bam" \
		M="${output_bam_dir}/${sample_ID}.re.sorted.du.matrix.txt" \
		REMOVE_DUPLICATES=false \
        CREATE_INDEX=true 2>> MARK_DUP.log

    # step3. 校正base品質
	gatk BaseRecalibrator \
		-R ${ref_fa} \
		-I "${output_bam_dir}/${sample_ID}.re.sorted.du.bam" \
		-O "${output_bam_dir}/${sample_ID}.recalibration_report.txt" \
		--known-sites ${dbsnp}

    # step4. 套用校正
    gatk ApplyBQSR \
		-R ${ref_fa} \
		-I "${output_bam_dir}/${sample_ID}.re.sorted.du.bam" \
		-O "${output_bam_dir}/${sample_ID}.re.sorted.du.bqsr.bam" \
		--bqsr-recal-file "${output_bam_dir}/${sample_ID}.recalibration_report.txt"
}

# config 
env_name="wes" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}
threads=20
picard="/home/${USER}/miniconda3/envs/${env_name}/share/picard-*/picard.jar" # 抓取使用者的環境資料夾


# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"

input_bam_dir="./output_bam_biobambam"
output_bam_dir="./bqsr_bam"
mkdir -p ./bqsr_bam
mkdir -p ./bqsr_bam/log

# col 6 sample_ID TCGA-4G-AAZR sort 並去除重複
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# for sample sample_IDs

for sample_ID in "${sample_IDs[@]}"
do

    echo "${sample_ID}"
	tumor="${sample_ID}-01A"
    normal="${sample_ID}-11A"
	
    bqsr ${tumor} >> "./bqsr_bam/log/${tumor}" 2>&1
    bqsr ${normal} >> "./bqsr_bam/log/${normal}" 2>&1
done