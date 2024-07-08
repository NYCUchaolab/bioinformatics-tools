#!/bin/bash
env_name="RNA_seq"
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh"
conda activate ${env_name}

cd /home/data/dataset/sra_samples/
gtf_file="/home/data/dataset/STAR_index/gencode.v22.annotation.gtf"
threads=20

samples=("SRR19117803" "SRR19117819" "SRR19117865")

for sample in "${samples[@]}"
do
    stringtie \
        -e -p 16 -G ${gtf_file} -o ./stringTie/${sample}_quant.gtf \
        ./HISAT2/${sample}.sorted.bam 1>./stringTie/${sample}_strg_quant.log 2>&1
        
    echo "${sample} ./stringTie/${sample}_quant.gtf" > "./stringTie/merge_list.txt"#³o¨B¬µ¤F
done


python prepDE.py3 -i ./stringTie/merge_list.txt -g .stringTie/output.count_matrix.csv -t .stringTie/output.transcript_count_matrix.csv