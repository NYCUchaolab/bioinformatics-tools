#!/bin/bash

# env
env_name="rna" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

# config
cd /home/data/dataset/sra_samples
mkdir -p ./RSEM_result
RSEM_index_dir="/home/data/dataset/RSEM_index"
STAR_output_dir="./2pass_0408"
RSEM_output_dir="./RSEM_result_0408"
threads=25

sample_list=("SRR19117803" "SRR19117819" "SRR19117865")
for sample in "${sample_list[@]}"
do
    rsem-calculate-expression \
        "${STAR_output_dir}/${sample}_Aligned.toTranscriptome.out.bam" \
        ${RSEM_index_dir} \
        "${RSEM_output_dir}/${sample}_" \
        -p $threads \
        --no-bam-output \
        --paired-end \
        --bam >> ./RSEM_result/log.txt
done

