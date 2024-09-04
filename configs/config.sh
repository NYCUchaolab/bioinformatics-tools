project_dir="/home/blue/code/test_variant_calling"
sample_sheet="${project_dir}/2sample_sheet.tsv"
GDC_token="${project_dir}/GDC_token.txt"
tools_dir="${project_dir}/bioinformatics-tools"
threads=8

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
germline="${db_dir}/af-only-gnomad.hg38.vcf.gz"
pon="${db_dir}/1000g_pon.hg38.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"
small_exac="${db_dir}/small_exac_common_3.hg38.vcf.gz"

# VEP setting
assembly=GRCh38

source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 