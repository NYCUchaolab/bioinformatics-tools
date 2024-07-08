#!/bin/bash
env_name="RNA_seq"
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh"
conda activate ${env_name}

cd /home/data/dataset/sra_samples
gtf_file="/home/data/dataset/STAR_index/gencode.v22.annotation.gtf"
threads=20

samples=("SRR19117803" "SRR19117819" "SRR19117865")

featureCounts \
    -T 20 -f -t gene \
    -g gene_id -a ${gtf_file} \
    -o ./featureCounts/all_sample.count \
    -p \
    ./HISAT2/${samples[0]}.sorted.bam \
    ./HISAT2/${samples[1]}.sorted.bam \
    ./HISAT2/${samples[2]}.sorted.bam