#!/bin/bash
download_GDC_files(){
    GDC_IDs=${1}
    bioinformatics-tools/gdc-client download ${GDC_IDs[*]} \
        -n ${threads} \
        -d "${download_dir}/tmp" \
        -t ${GDC_token} \
        --log-file "${download_dir}/log/download.log" \
        --latest 
}

# conda 
env_name="wes-preprocessing" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

# config 
source "../config.sh"
cd ${project_dir}

# mkdir
download_dir="${project_dir}/bam"
mkdir -p "${download_dir}"
mkdir -p "${download_dir}/tmp"
mkdir -p "${download_dir}/log"

# load sample sheet 
GDC_IDs=($(tail -n +2 ${sample_sheet} | cut -f 1 | sort -u))
download_GDC_files ${GDC_IDs}

# move and rename
tail -n +2 $sample_sheet | while IFS=$'\t' read -r -a sample
do
    GDC_ID=${sample[0]}
    file_name=${sample[1]}
    new_sample_ID=${sample[8]}

    bam_path=${download_dir}/tmp/${GDC_ID}/${file_name}
    mv ${bam_path} ${download_dir}/${new_sample_ID}.bam
	mv ${bam_path/.bam/.bai} ${download_dir}/${new_sample_ID}.bam.bai
done

# delete files
#rm -r "${download_dir}/tmp"
python "${tools_dir}/others/send_msg.py" "download_GDC_files done"