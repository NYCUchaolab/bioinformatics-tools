#!/bin/bash

env_name="rna" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

cd /home/data/dataset/sra_samples
star_1pass_dir="./1pass_0408"
star_2pass_dir="./2pass_0408"
RSEM_result_dir="./RSEM_result_0408"

# ref
STAR_index_dir="/home/data/dataset/STAR_index"
kallisto_index_dir="/home/data/dataset/kallisto_index"
RSEM_index_dir="/home/data/dataset/RSEM_index"
threads=25

mkdir -p ${star_1pass_dir}
mkdir -p ${star_2pass_dir}
mkdir -p ${RSEM_index_dir}
mkdir -p ${RSEM_result_dir}

# create star 1pass index (可直接抓取做好的)
# check 是否存在
sample_list=("SRR19117803" "SRR19117819" "SRR19117865")

# 1-pass index
#STAR \
#    --runMode genomeGenerate \
#   --genomeDir ${star_1pass_dir} \
#    --sjdbGTFfile ${STAR_index_dir}/gencode.v22.annotation.gtf \
#    --genomeFastaFiles ${STAR_index_dir}/GRCh38.d1.vd1.fa \
#    --runThreadN $threads \
#    --outFileNamePrefix ${star_1pass_dir}/ > "${star_1pass_dir}/log.txt"

# 1-pass mapping
for sample in "${sample_list[@]}"
do
    STAR \
        --genomeDir ${STAR_index_dir} \
        --readFilesIn ./fastq/${sample}_1.fastq ./fastq/${sample}_2.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts TranscriptomeSAM \
        --runThreadN $threads \
        --outFileNamePrefix ${star_1pass_dir}/${sample}_ >> "${star_1pass_dir}/log.txt"
done

# 1-pass output path
output_list_1pass=() 
for sample in "${sample_list[@]}"; do
    output_list_1pass+=( "${star_1pass_dir}/${sample}_SJ.out.tab" )
done

# 2-pass index
STAR \
    --runMode genomeGenerate \
    --genomeDir ${star_2pass_dir} \
    --genomeFastaFiles ${STAR_index_dir}/GRCh38.d1.vd1.fa \
    --runThreadN $threads \
    --sjdbGTFfile ${STAR_index_dir}/gencode.v22.annotation.gtf \
    --sjdbFileChrStartEnd "${output_list_1pass[@]}" \
    --outFileNamePrefix ${star_2pass_dir}/ > "${star_2pass_dir}/log.txt"

# 2-pass mapping
for sample in "${sample_list[@]}"
do  
    STAR \
        --genomeDir ${star_2pass_dir} \
        --readFilesIn ./fastq/${sample}_1.fastq ./fastq/${sample}_2.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts TranscriptomeSAM \
        --runThreadN $threads \
        --outFileNamePrefix ${star_2pass_dir}/${sample}_ >> "${star_2pass_dir}/log.txt"
done

for sample in "${sample_list[@]}"
do
    rsem-calculate-expression \
        ${star_2pass_dir}/${sample}_Aligned.toTranscriptome.out.bam \
        ${RSEM_index_dir} \
        ${RSEM_result_dir}/${sample}_ \
        -p $threads --no-bam-output --paired-end --bam

    kallisto quant \
        -i ${kallisto_index_dir}/kallisto_index.idx \
        -o kallisto/${sample} \
        -t $threads \
        ./fastq/${sample}_1.fastq \
        ./fastq/${sample}_2.fastq
done
