#!/bin/bash
env_name="RNA_seq"
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh"
conda activate ${env_name}

cd /home/data/dataset/sra_samples/

threads=20
#index

salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i athal_index_salmon

#run
samples=("SRR19117803" "SRR19117819" "SRR19117865")

for sample in "${samples[@]}"
do
    salmon quant -i ./kallisto_salmon/athal_index_salmon -l A \
    		-1 ./fastq/${sample}_1.fastq \
    		-2 ./fastq/${sample}_2.fastq \
    		-p ${threads} \
        -o ./kallisto_salmon/${sample}_quant_
done

