#!/bin/bash
run_varscan(){

    sample_ID=${1}
    tumor="${case_ID}-01A"
    normal="${case_ID}-11A"

    # Step 1: Mpileup; Samtools
    samtools mpileup \
    -f ${ref_fa} \
    -q 1 \
    -B  \
    ./sorted_du_bqsr_bam_biobambam/${normal}.sorted.du.bqsr.bam \
    ./sorted_du_bqsr_bam_biobambam/${tumor}.sorted.du.bqsr.bam \
    > ./variants_calling/varscan/${sample_ID}_intermediate_mpileup.pileup
    echo "step1 samtools mpileup ${sample_ID} done"
    
    # Step 2: Varscan Somatic; Varscan.v2
    java -jar ${varscan} somatic \
    "./variants_calling/varscan/${sample_ID}_intermediate_mpileup.pileup" \
    "./variants_calling/varscan/${sample_ID}.varscan" \
    --mpileup 1 \
    --min-var-freq 0.1 \
    --output-vcf 
    echo "step2 somatic ${sample_ID} done"
    
    # Step 3: Varscan ProcessSomatic; Varscan.v2
    java -jar ${varscan} processSomatic \
    "./variants_calling/varscan/${sample_ID}.varscan.snp.vcf" \
    --min-tumor-freq 0.10 \
    --max-normal-freq 0.05 \
    --p-value 0.07
    
    java -jar ${varscan} processSomatic \
    "./variants_calling/varscan/${sample_ID}.varscan.indel.vcf" \
    --min-tumor-freq 0.1 \
    --max-normal-freq 0.05 \
    --p-value 0.07
    echo "step3 processSomatic ${sample_ID} done"

    # Step 4: del tmp
    rm "./variants_calling/varscan/${sample_ID}_intermediate_mpileup.pileup"
    echo "step4 clean ${sample_ID} tmp file done"

}

# conda 
env_name="wes" # 環境名稱請確認
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda activate ${env_name}

# ref
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa" # 修改ref位置
varscan="/home/${USER}/miniconda3/envs/${env_name}/share/varscan-*/VarScan.jar" # 抓取使用者的環境資料夾

# config 
project_dir="/home/data/dataset/CHOL_10sample" # 專案資料夾
cd ${project_dir}

# mkdir
mkdir -p ./variants_calling/varscan && chmod a+rw ./variants_calling/varscan

# load sample sheet 
sample_sheet="./s5sheet.tsv"
sample_IDs=($(tail -n +2 ${sample_sheet} | cut -f 6 | sort -u))

# for sampleID in sample_IDs
for sample_ID in "${sample_IDs[@]}"
do
    run_varscan ${sample_ID} &
done


