#!/bin/bash

# conda 
env_name="wes" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}
# cnvkit="$HOME/miniconda3/envs/wes/bin/cnvkit.py"

# ref
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa" # 修改ref位置

# config 
threads=20
project_dir="/home/data/dataset/CHOL_10sample" # 專案資料夾
cd ${project_dir}
access="./cnv_ref/access-10kb.hg38.bed"
refFlat="./cnv_ref/refFlat.txt"
#illumina_bed="./cnv_ref/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"
illumina_bed="./cnv_ref/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.bed"
input_bam_dir="sorted_du_bqsr_bam_biobambam"
cnv_output_dir="./variants_calling/cnv"

# mkdir
mkdir -p variants_calling
mkdir -p ${cnv_output_dir}
mkdir -p ${cnv_output_dir}/log
cd ${project_dir}

cnvkit.py batch ${project_dir}/${input_bam_dir}/*01A*.bam \
        --normal ${project_dir}/${input_bam_dir}/*11A*.bam \
        --targets ${illumina_bed} --annotate ${refFlat} \
        --fasta ${ref_fa} --access ${access} -p ${threads}\
        --output-reference ${cnv_output_dir}/cnv_ref.cnn --output-dir ${cnv_output_dir} \
        --diagram --scatter  > ${cnv_output_dir}/log/cnvkit.log 2>&1