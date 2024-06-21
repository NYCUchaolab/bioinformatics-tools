#!/bin/bash

run_cnvkit(){
        bam_input_dir=${1}
        cnv_output_dir=${2}
        # call copy number variants
        ${cnvkit} batch ${project_dir}/${bam_input_dir}/*T*.bam \
                --normal ${project_dir}/${bam_input_dir}/*B*.bam \
                --targets ${illumina_bed} --annotate ${refFlat} \
                --fasta ${ref_fa} --access ${access} -p ${threads}\
                --output-reference ${cnv_output_dir}/cnv_ref.cnn --output-dir ${cnv_output_dir} \
                --diagram --scatter  > ${cnv_output_dir}/log/cnvkit.log 2>&1
        # export to gistic format
        ${cnvkit} export seg $cnv_output_dir/*.cns \
                -o ${cnv_output_dir}/gistic.segments >> ${cnv_output_dir}/log/cnvkit.log 2>&1
}

run_gistic2(){
        gistic_output_dir={$1}
        # run gistic
        $gistic -b ${gistic_output_dir} \
                -seg ${cnv_output_dir}/gistic.segments\
                -refgene ${gistic_ref} \
                -genegistic 1 \
                -smallmem 1 \
                -broad 1 \
                -brlen 0.5 \
                -conf 0.90 -armpeel 1 -savegene 1 -gcm "extreme"
}

# conda 
env_name="wes-cnv" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}
cnvkit="$HOME/miniconda3/envs/${env_name}/bin/cnvkit.py"

# ref
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa" # 修改ref位置
gistic_dir="/home/data/data_dingyangliu/ntu-tki/GISTIC2"
gistic_ref="${gistic_dir}/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"
gistic="${gistic_dir}/gistic2"
#export PATH=${gistic_dir}:$PATH
export MCR_ROOT=/home/data/data_dingyangliu/ntu-tki/GISTIC2/MATLAB_Compiler_Runtime
export LD_LIBRARY_PATH=$MCR_ROOT/v92/runtime/glnxa64:$MCR_ROOT/v92/bin/glnxa64:$MCR_ROOT/v92/sys/os/glnxa64:$LD_LIBRARY_PATH

# config 
project_dir="/home/data/data_dingyangliu/ntu-tki" # 專案資料夾
cd ${project_dir}
access="access-10kb.hg38.bed"
refFlat="refFlat_hg38.txt"
illumina_bed="hg38_Twist_ILMN_Exome_2.5_Panel_Combined_Mito.bed"

bam_input_dir="bam"
cnv_output_dir="./variants_calling/cnv_TB"
gistic_output_dir="./variants_calling/gistic_test_TB"

# mkdir
mkdir -p variants_calling
mkdir -p ${cnv_output_dir}
mkdir -p ${cnv_output_dir}/log
mkdir -p ${gistic_output_dir}
cd ${project_dir}

# call functions 
# run_cnvkit ${bam_input_dir} ${cnv_output_dir}
run_gistic2 ${gistic_output_dir}