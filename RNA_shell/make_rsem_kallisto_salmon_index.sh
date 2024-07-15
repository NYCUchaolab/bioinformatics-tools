#!/bin/bash

# env
env_name="rna" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

# config
STAR_index_dir="/home/data/dataset/STAR_index"
RSEM_index_dir="/home/data/dataset/rsem_index"
kallisto_index_dir="/home/data/dataset/kallisto_index"
threads=25

# make rsem-index
rsem-prepare-reference \
    ${STAR_index_dir}/GRCh38.d1.vd1.fa \
    ${RSEM_index_dir}/ \
    --gtf ${STAR_index_dir}/gencode.v22.annotation.gtf \
    -p {threads}

# download cdna.fa file
rm "${kallisto_index_dir}/Homo_sapiens.GRCh38.cdna.all.fa.gz"
wget -P ${kallisto_index_dir} https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

cd ${kallisto_index_dir}

# make kallisto-index
kallisto index \
    -i "${kallisto_index_dir}/kallisto_index.idx" \
    "${kallisto_index_dir}/Homo_sapiens.GRCh38.cdna.all.fa.gz"

# make salmon-index
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i athal_index_salmon