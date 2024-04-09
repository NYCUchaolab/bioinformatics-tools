input_vcf="/home/data/dataset/CHOL_10sample/variants_calling/mutect/TCGA-4G-AAZF.vcf"
output_vcf="/home/data/dataset/CHOL_10sample/variants_calling/mutect/TCGA-4G-AAZF_splited.vcf"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"

bcftools norm \
    -m-any \
    --check-ref \
    -xw \
    -f ${ref_fa} \
    ${input_vcf} \
    -o ${output_vcf}
