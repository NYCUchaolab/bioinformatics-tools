#!/bin/bash
env_name="rna" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

cd /home/data/dataset/sra_samples
kallisto_index_dir="/home/data/dataset/kallisto_index"
mkdir -p kallisto

threads=25
samples=("SRR19117803" "SRR19117819" "SRR19117865")

for sample in "${samples[@]}"
do
    kallisto quant \
        -i ${kallisto_index_dir}/kallisto_index.idx \
        -o kallisto/${sample} \
        -t $threads \
        ./fastq/${sample}_1.fastq \
        ./fastq/${sample}_2.fastq
done
