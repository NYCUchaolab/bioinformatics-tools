#!/bin/bash
run_samtools_sorting(){
    #samtools sort \
        #    -n -@ ${threads} \
        #    -o ${bam_dir}/${sample_ID}.sorted.bam \
        #    ${bam_dir}/${sample_ID}.bam > ${bam_dir}/log/${sample_ID}.sort.log
}

run_bam2fastq() {
    sample_ID=${1}
    # 注意高記憶體用量
	
	bedtools bamtofastq \
        -i ${bam_dir}/${sample_ID}.sorted.bam  \
		-fq ${fastq_dir}/${sample_ID}_R1.fq \
		-fq2 ${fastq_dir}/${sample_ID}_R2.fq > ${fastq_dir}/log/${sample_ID}.bam2fastq.log
    # run fastQC 
    fastqc -o "${qc_dir}/fastqc" "${fastq_dir}/${sample_ID}_R1.fq"
    fastqc -o "${qc_dir}/fastqc" "${fastq_dir}/${sample_ID}_R2.fq"
    python rna_shell/send_msg.py "${sample_ID} done"
}

# config 
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate wes


# project config
source "../config.sh"
cd ${project_dir}
bam_dir="${project_dir}/bam"
fastq_dir="${project_dir}/fastq"
qc_dir="${project_dir}/qc"
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 9 | sort -u))

# package config and mkdir
mkdir -p ${fastq_dir}
mkdir -p ${fastq_dir}/log
mkdir -p ${bam_dir}/log
mkdir -p ${qc_dir}
mkdir -p ${qc_dir}/fastqc
mkdir -p ${qc_dir}/multiqc

pids=()
# run bam2fastq
for sample_ID in "${sample_IDs[@]}"
do  
    echo ${sample_ID}
    # sort bam and bam2fastq
    run_bam2fastq "${sample_ID}" &
    pids+=($!)
done

# run multiQC
for pid in "${pids[@]}"
do
    wait "$pid"
done
multiqc -o "${qc_dir}/multiqc" "${qc_dir}/fastqc"
python send_msg.py "bam2fastq_and_QC done"