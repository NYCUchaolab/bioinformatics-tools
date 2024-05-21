project_dir="/home/data/dataset/CHOL_2sample_new"
sample_sheet="${project_dir}/gdc_sample_sheet_test.2024-02-26.tsv"
GDC_token="/home/data/tools/gdc-user-token.2024-05-15T09_11_50.917Z.txt"
tools_dir="${project_dir}/bioinformatics-tools"
threads=20

# ref config
db_dir="/home/data/dataset/gatk_file"
dbsnp="${db_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
germline="${db_dir}/af-only-gnomad.hg38.vcf.gz"
pon="${db_dir}/1000g_pon.hg38.vcf.gz"
ref_fa="/home/data/dataset/BWA_index/GRCh38.d1.vd1.fa"
small_exac="${db_dir}/small_exac_common_3.hg38.vcf.gz"

source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 