#!/bin/bash
env_name="RNA_seq"
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh"
conda activate ${env_name}

cd /home/data/dataset/sra_samples/

threads=20

samples=("SRR19117803" "SRR19117819" "SRR19117865")

for sample in "${samples[@]}"
do
    hisat2 -x /home/data/dataset/sra_samples/HISAT2/ref/grch38/genome \
        -p 20 \
        -1 ./fastq/${sample}_1.fastq \
        -2 ./fastq/${sample}_2.fastq \
        -S ./HISAT2/${sample}.sam > ./HISAT2/${sample}.log
    
    samtools sort -@ 20 -o ./HISAT2/${sample}.sorted.bam ./HISAT2/${sample}.sam
    
    samtools index ./HISAT2/${sample}.sorted.bam
    
    samtools flagstat ./HISAT2/${sample}.sorted.bam > ./HISAT2/${sample}.flagstat.log
done