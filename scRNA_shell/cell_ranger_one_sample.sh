#!/bin/bash

# add spaceranger to PATH
export PATH=/home/data/tools/cellranger-8.0.0:$PATH

# configs
output_dir="rna_test"
core=15
ram=64
fastq_dir="new_fastq"
sample_ID="SRR6334436"
cellrange_ref="/home/data/dataset/refdata-gex-GRCh38-2020-A"

# make dir
mkdir -p ${output_dir}

# cellranger
cellranger count --id=${output_dir}\
    --localcores=${core} \ 
    --localmem=${ram} \
    --sample=${sample_ID} \
    --transcriptome=${cellrange_ref} \
    --fastqs=${fastq_dir} \
    --create-bam=true \
    > "${output_dir}/log_file.log" 2>&1
