#!/bin/bash

run_muse(){

    sample_ID=${1}
    
    tumor_file="${bam_dir}/${sample_ID}-01A.sorted.du.bqsr.bam"
    normal_file="${bam_dir}/${sample_ID}-11A.sorted.du.bqsr.bam"
    
    MuSE call -f ${ref_fa} -O ${muse_dir}/${sample_ID} ${tumor_file} ${normal_file} > "${muse_dir}/log/${sample_ID}.log" 2>&1 &
    MuSE_pid=$! # 獲取最後一個執行的pid
    wait "${MuSE_pid}" # 等待muse 結束

    MuSE_exit_status=$?
    if [ "${MuSE_exit_status}" -eq 0 ]; then
        echo "${sample_ID} MuSE call completed successfully" >> "${muse_dir}/muse.log"
    else
        echo "${sample_ID} MuSE call failed with exit status ${MuSE_exit_status}" >> "${muse_dir}/muse.log"
    fi

    MuSE sump -I ${muse_dir}/${sample_ID}.MuSE.txt -E -D ${dbsnp} -O ${muse_dir}/${sample_ID}.MuSE.vcf  >> "${muse_dir}/log/${sample_ID}.log" 2>&1 &
}

# config 
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate wes

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet.2024-02-26.tsv"
bam_dir="${project_dir}/sorted_du_bqsr_bam_biobambam"

# package config and mkdir
muse_dir="variants_calling/MuSE"
mkdir -p variants_calling && chmod a+rw variants_calling
mkdir -p variants_calling/MuSE && chmod a+rw variants_calling/MuSE
mkdir -p variants_calling/MuSE/log && chmod a+rw variants_calling/MuSE/log

# col 6 sample_ID TCGA-4G-AAZR sort並去除重複
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# cp files
for sample_ID in "${sample_IDs[@]}"
do
    #echo "${sample_ID}"
    normal="${sample_ID}-01A"
    tumor="${sample_ID}-11A"
    cp -n "${bam_dir}/${normal}.sorted.du.bqsr.bai" "${bam_dir}/${normal}.sorted.du.bqsr.bam.bai"
    wait
    cp -n "${bam_dir}/${tumor}.sorted.du.bqsr.bai" "${bam_dir}/${tumor}.sorted.du.bqsr.bam.bai"
done

# run muse
for sample_ID in "${sample_IDs[@]}"
do  
    echo ${sample_ID}
    run_muse "${sample_ID}" &
done
