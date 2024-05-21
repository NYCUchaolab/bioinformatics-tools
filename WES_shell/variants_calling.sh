#!/bin/bash

run_mutect(){
    case_ID=${1}
    tumor_file="${bam_dir}/${case_ID}-T.sorted.du.bqsr.bam"
    normal_file="${bam_dir}/${case_ID}-N.sorted.du.bqsr.bam"
    
    gatk --java-options "-Xmx4G" Mutect2 \
    -R ${ref_fa} \
    -I ${tumor_file} \
    -I ${normal_file} \
    -tumor "${case_ID}-T" \
    -normal "${case_ID}-N" \
    --germline-resource ${germline} \
    --native-pair-hmm-threads ${threads} \
    -pon ${pon} \
    -O "${mutect_dir}/${case_ID}.vcf" > "${mutect_dir}/log/${case_ID}.log" 2>&1
    
    gatk --java-options "-Xmx10G" GetPileupSummaries \
    -I ${tumor_file} \
    -V ${small_exac} \
    -L ${small_exac}  \
    -O "${mutect_dir}/${case_ID}-T.pileups.table" >> "${mutect_dir}/log/${case_ID}.log" 2>&1

    gatk --java-options "-Xmx10G" GetPileupSummaries \
    -I ${normal_file} \
    -V ${small_exac} \
    -L ${small_exac}  \
    -O "${mutect_dir}/${case_ID}-N.pileups.table" >> "${mutect_dir}/log/${case_ID}.log" 2>&1

    gatk --java-options "-Xmx10G" CalculateContamination \
    -I "${mutect_dir}/${case_ID}-T.pileups.table" \
    -matched "${mutect_dir}/${case_ID}-N.pileups.table" \
    -O "${mutect_dir}/${case_ID}.contamination.table" \
    -tumor-segmentation "${mutect_dir}/${case_ID}.segments.table" >> "${mutect_dir}/log/${case_ID}.log" 2>&1

    gatk --java-options "-Xmx10G" FilterMutectCalls \
    -V "${mutect_dir}/${case_ID}.vcf" \
    --tumor-segmentation "${mutect_dir}/${case_ID}.segments.table" \
    --contamination-table "${mutect_dir}/${case_ID}.contamination.table" \
    -O "${mutect_dir}/${case_ID}.filtered.vcf" \
    -R ${ref_fa} >> "${mutect_dir}/log/${case_ID}.log" 2>&1

    bcftools norm \
    -m-any \
    -cx\
    --check-ref \
    -xw \
    -f ${ref_fa} \
    "${mutect_dir}/${case_ID}.filtered.vcf" \
    -o "${mutect_dir}/${case_ID}.splited.vcf" >> "${mutect_dir}/log/${case_ID}.log" 2>&1
}


run_muse(){
    case_ID=${1}
    tumor_file="${bam_dir}/${case_ID}-T.sorted.du.bqsr.bam"
    normal_file="${bam_dir}/${case_ID}-N.sorted.du.bqsr.bam"
    
    MuSE call -f ${ref_fa} -O ${muse_dir}/${case_ID} ${tumor_file} ${normal_file} > "${muse_dir}/log/${case_ID}.MuSE.log" 2>&1 
    MuSE sump -I ${muse_dir}/${case_ID}.MuSE.txt -E -D ${dbsnp} -O ${muse_dir}/${case_ID}.MuSE.vcf  >> "${muse_dir}/log/${case_ID}.MuSE.log" 2>&1
    }

run_varscan(){

    case_ID=${1}
    tumor="${case_ID}-T"
    normal="${case_ID}-N"

    # Step 1: Mpileup; Samtools
    samtools mpileup \
    -f ${ref_fa} \
    -q 1 \
    -B \
    ${bam_dir}/${normal}.sorted.du.bqsr.bam \
    ${bam_dir}/${tumor}.sorted.du.bqsr.bam \
    > ${varscan_dir}/${case_ID}_intermediate_mpileup.pileup
    echo "step1 samtools mpileup ${case_ID} done"
    
    # Step 2: Varscan Somatic; Varscan.v2
    java -jar ${varscan} somatic \
    "${varscan_dir}/${case_ID}_intermediate_mpileup.pileup" \
    "${varscan_dir}/${case_ID}.varscan" \
    --mpileup 1 \
    --min-var-freq 0.1 \
    --output-vcf 
    echo "step2 somatic ${case_ID} done"
    
    # Step 3: Varscan ProcessSomatic; Varscan.v2
    java -jar ${varscan} processSomatic \
    "${varscan_dir}/${case_ID}.varscan.snp.vcf" \
    --min-tumor-freq 0.10 \
    --max-normal-freq 0.05 \
    --p-value 0.05
    
    java -jar ${varscan} processSomatic \
    "${varscan_dir}/${case_ID}.varscan.indel.vcf" \
    --min-tumor-freq 0.10 \
    --max-normal-freq 0.05 \
    --p-value 0.05
    echo "step3 processSomatic ${case_ID} done"

    # Step 4: del tmp
    rm "${varscan_dir}/${case_ID}_intermediate_mpileup.pileup"
    echo "step4 clean ${case_ID} tmp file done"

}

generate_pindel_config(){
    case_ID=${1}
    echo "${bam_dir}/${case_ID}-N.sorted.du.bqsr.bam 250 ${case_ID}_N" > ${config_file}
    echo "${bam_dir}/${case_ID}-T.sorted.du.bqsr.bam 250 ${case_ID}_T" >> ${config_file}    
    
    echo "indel.filter.input = ${pindel_dir}/${case_ID}.all.head" > ${filter_config_file}
    echo "indel.filter.vaf = 0.05" >> ${filter_config_file}
    echo "indel.filter.cov = 20" >> ${filter_config_file}
    echo "indel.filter.hom = 6" >> ${filter_config_file}
    echo "indel.filter.pindel2vcf = pindel2vcf" >> ${filter_config_file}
    echo "indel.filter.reference = ${ref_fa}" >> ${filter_config_file}
    echo "indel.filter.referencename = $(basename "${ref_fa}")" >> ${filter_config_file}
    echo "indel.filter.referencedate = $(date +%Y%m%d)" >> ${filter_config_file}
    echo "indel.filter.output = ${pindel_dir}/${case_ID}.indel.filtered.vcf" >> ${filter_config_file}
}


run_pindel(){
    case_ID=${1}
    config_file="${pindel_dir}/${case_ID}.config.txt"
    filter_config_file="${pindel_dir}/${case_ID}.indel.filter.config"
    generate_pindel_config ${case_ID}

    # step 1: pindel variant calling
    pindel -f ${ref_fa} -i ${config_file} -c ALL -T ${threads} -o ${pindel_dir}/${case_ID}.Pindel

    # step 2: extract indel summary lines
    grep ChrID ${pindel_dir}/${case_ID}.Pindel_D > ${pindel_dir}/${case_ID}.D.head
    grep ChrID ${pindel_dir}/${case_ID}.Pindel_SI > ${pindel_dir}/${case_ID}.SI.head
    cat ${pindel_dir}/${case_ID}.D.head ${pindel_dir}/${case_ID}.SI.head > ${pindel_dir}/${case_ID}.all.head

    #step 3: run the provided perl script with an updated configuration file
    #cd ${pindel_dir}
    
    perl "${tools_dir}/WES_shell/somatic_indelfilter.pl" "${pindel_dir}/${case_ID}.indel.filter.config"
    
    # ${case_ID} 
    rm ${config_file}
    rm ${filter_config_file}
    
}
