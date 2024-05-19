#!/bin/bash
run_sort(){
    sample_ID=${1}
    rm ${bam_dir}/${sample_ID}.sorted.bam.tmp*
    samtools sort \
            -n -@ ${threads} \
            -o ${bam_dir}/${sample_ID}.sorted.bam \
            ${bam_dir}/${sample_ID}.bam > ${bam_dir}/log/${sample_ID}.sort.log 2>&1
}

run_bam2fastq() {
	bedtools bamtofastq \
        -i ${bam_dir}/${sample_ID}.sorted.bam  \
		-fq ${fastq_dir}/${sample_ID}_R1.fq \
		-fq2 ${fastq_dir}/${sample_ID}_R2.fq > ${fastq_dir}/log/${sample_ID}.bam2fastq.log 2>&1
    # run fastQC 
    fastqc -o "${qc_dir}/fastqc" "${fastq_dir}/${sample_ID}_R1.fq" > ${fastq_dir}/${sample_ID}_R1.fastqc.log 2>&1
    fastqc -o "${qc_dir}/fastqc" "${fastq_dir}/${sample_ID}_R2.fq" > ${fastq_dir}/${sample_ID}_R2.fastqc.log 2>&1

    python "${tools_dir}/others/send_msg.py" "${sample_ID} bamtofastq done"
}


# conda 
env_name="wes-preprocessing" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

# config 
source "../config.sh"

# project config
cd ${project_dir}
bam_dir="${project_dir}/bam"
fastq_dir="${project_dir}/fastq"
qc_dir="${project_dir}/qc"

# package config and mkdir
mkdir -p ${fastq_dir}
mkdir -p ${fastq_dir}/log
mkdir -p ${bam_dir}/log
mkdir -p ${qc_dir}
mkdir -p ${qc_dir}/fastqc
mkdir -p ${qc_dir}/multiqc

# load sample sheet 
sample_IDs=($(tail -n +2 ${sample_sheet} | cut -f 9 | sort -u))
for sample_ID in "${sample_IDs[@]}"
    do  
        run_sort ${sample_ID}
    done

pids=()
for sample_ID in "${sample_IDs[@]}"
    do  
        echo ${sample_ID}
        # sort bam and bam2fastq
        run_bam2fastq "${sample_ID}" &
        pids+=($!)
    done

# wait all file done
for pid in "${pids[@]}"
    do
        wait "$pid"
    done

multiqc --force -o "${qc_dir}/multiqc" "${qc_dir}/fastqc" 
python "${tools_dir}/others/send_msg.py" "bam2fastq_and_QC done"