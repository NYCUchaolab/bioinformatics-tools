#!/bin/bash

run_mutect(){
    case_ID=${1}
    tumor_file="${bam_dir}/${case_ID}-01A.re.sorted.du.bqsr.bam" # 注意後綴
    normal_file="${bam_dir}/${case_ID}-11A.re.sorted.du.bqsr.bam" # 注意後綴
    
    gatk Mutect2 \
    -R ${ref_fa} \
    -I ${tumor_file} \
    -I ${normal_file} \
    -normal "${case_ID}-11A" \
    -tumor "${case_ID}-01A" \
    --germline-resource ${germline} \
    -pon ${pon} \
    -O "${mutect_dir}/${case_ID}.vcf" > "${mutect_dir}/log/${case_ID}.log" 2>&1
    
    gatk GetPileupSummaries \
    -I ${tumor_file} \
    -V ${small_exac} \
    -L ${small_exac}  \
    -O "${mutect_dir}/${case_ID}-01A.pileups.table" >> "${mutect_dir}/log/${case_ID}.log" 2>&1

    gatk GetPileupSummaries \
    -I ${normal_file} \
    -V ${small_exac} \
    -L ${small_exac}  \
    -O "${mutect_dir}/${case_ID}-11A.pileups.table"

    gatk CalculateContamination \
    -I "${mutect_dir}/${case_ID}-01A.pileups.table" \
    -matched "${mutect_dir}/${case_ID}-11A.pileups.table" \
    -O "${mutect_dir}/${case_ID}.contamination.table" \
    -tumor-segmentation "${mutect_dir}/${case_ID}.segments.table" >> "${mutect_dir}/log/${case_ID}.log" 2>&1

    gatk FilterMutectCalls \
    -V "${mutect_dir}/${case_ID}.vcf" \
    --tumor-segmentation "${mutect_dir}/${case_ID}.segments.table" \
    --contamination-table "${mutect_dir}/${case_ID}.contamination.table" \
    -O "${mutect_dir}/${case_ID}.filtered.vcf" \
    -R ${ref_fa} >> "${mutect_dir}/log/${case_ID}.log" 2>&1

}


# config 
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate wes

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
germline="${db_dir}/af-only-gnomad.hg38.vcf.gz"
pon="${db_dir}/1000g_pon.hg38.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

# project config
project_dir="/home/data/dataset/CHOL_10sample"
cd ${project_dir}
samplesheet_path="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"
bam_dir="${project_dir}/bqsr_bam"

# package config and mkdir
mutect_dir="variants_calling/mutect"
mkdir -p variants_calling && chmod a+rw variants_calling
mkdir -p variants_calling/mutect && chmod a+rw variants_calling/mutect
mkdir -p variants_calling/mutect/log && chmod a+rw variants_calling/mutect/log

# col 6 sample_ID TCGA-4G-AAZR sort並去除重複
sample_IDs=($(tail -n +2 ${samplesheet_path} | cut -f 6 | sort -u))

# run muse
for sample_ID in "${sample_IDs[@]}"
do  
    echo ${sample_ID}
    run_mutect "${sample_ID}"
done
